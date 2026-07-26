"""Ligand preparation: compound library in, docking-ready PDBQT out.

The pipeline per compound is:

1. **Parse** the record with RDKit and sanitise it.
2. **Protonate** at the target pH, and optionally enumerate tautomers.
3. **Embed and minimise** in 3D (ETKDG followed by MMFF94s).
4. **Convert** to PDBQT with Meeko, which assigns AutoDock atom types, computes
   partial charges and decides which bonds are rotatable.
5. **Describe** the compound so the eventual hit list is interpretable.

Steps 2 and 3 are both done by Scrubber in a single pass, which is also where the
embedding seed is set -- see :func:`prepare_one`.

A note on ring conformations
----------------------------

Vina samples rotatable bonds but never ring geometry: whatever ring conformation
is in the PDBQT is held rigid for the whole docking run. A sugar that arrives as a
strained boat stays a boat.

MMFF94s minimisation (the default) relaxes each molecule into its *nearest*
energy minimum, which repairs bad input geometry but will not cross a barrier to
find a better ring pucker. ``n_conformers > 1`` embeds several distinct starting
geometries and docks each, which is the only setting here that genuinely samples
ring space. Meeko's macrocycle handling -- on by default -- separately lets Vina
open and re-close large rings during the search.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from drydock.core.descriptors import Descriptors2D
from drydock.core.library import Record

# Number of active torsions Vina can handle before the search becomes
# unproductive. Not a hard engine limit, but past roughly this many the pose is
# unlikely to mean much, and it is worth telling the user rather than burning
# hours on it.
TORSION_WARN_THRESHOLD = 32


@dataclass(frozen=True, slots=True)
class PrepConfig:
    """Settings for a preparation run. Recorded in provenance verbatim."""

    ph: float = 7.4
    """Target pH for protonation state assignment."""

    optimize: bool = True
    """Embed in 3D and minimise with MMFF94s. Disable to keep input coordinates."""

    n_conformers: int = 1
    """Starting geometries per compound. >1 is the only real ring sampling."""

    skip_tautomers: bool = True
    """Enumerating tautomers multiplies the library; off by default."""

    macrocycles: bool = True
    """Let Vina open and re-close macrocyclic rings during the search."""

    seed: int = 0
    """Global seed; per-ligand seeds are derived from it deterministically."""

    max_torsions: int | None = None
    """Reject ligands with more rotatable bonds than this. None disables."""

    def to_dict(self) -> dict[str, Any]:
        from dataclasses import asdict

        return asdict(self)


@dataclass(frozen=True, slots=True)
class PreparedLigand:
    """A successfully prepared ligand, one manifest row."""

    ligand_id: str
    compound_id: str
    pdbqt_files: tuple[str, ...]
    torsions: int
    prepared_charge: int
    descriptors: Descriptors2D
    warnings: tuple[str, ...] = ()

    def to_row(self) -> dict[str, Any]:
        """Flatten into the manifest's column layout."""
        row: dict[str, Any] = {
            "ligand_id": self.ligand_id,
            "compound_id": self.compound_id,
            "n_conformers": len(self.pdbqt_files),
            "pdbqt": ";".join(self.pdbqt_files),
            "torsions": self.torsions,
            "prepared_charge": self.prepared_charge,
        }
        row.update(self.descriptors.to_dict())
        if self.warnings:
            row["warnings"] = "; ".join(self.warnings)
        return row


@dataclass(frozen=True, slots=True)
class PrepFailure:
    """A compound that could not be prepared, and why."""

    ligand_id: str
    compound_id: str
    stage: str
    error: str

    def to_row(self) -> dict[str, Any]:
        return {
            "ligand_id": self.ligand_id,
            "compound_id": self.compound_id,
            "stage": self.stage,
            "error": self.error,
        }


@dataclass
class _Tools:
    """Per-process RDKit/Meeko/Scrubber objects.

    Building a Scrub instance parses pKa and tautomer definition files, which is
    slow enough that doing it per ligand would dominate preparation. Workers
    build one on first use and keep it.
    """

    scrub: Any = None
    scrub_keep_coords: Any = None
    preparator: Any = None
    config: PrepConfig | None = field(default=None)


_TOOLS = _Tools()


def _get_tools(config: PrepConfig) -> _Tools:
    """Return this process's cached tools, rebuilding if the config changed."""
    if _TOOLS.config == config and _TOOLS.preparator is not None:
        return _TOOLS

    from meeko import MoleculePreparation
    from molscrub import Scrub
    from rdkit import RDLogger

    # RDKit is voluble about valence and stereochemistry issues that we handle by
    # recording a failure. Left enabled it interleaves with progress output from
    # every worker at once.
    RDLogger.DisableLog("rdApp.*")

    _TOOLS.scrub = Scrub(
        ph_low=config.ph,
        ph_high=None,
        skip_tautomers=config.skip_tautomers,
        skip_gen3d=not config.optimize,
        numconfs=config.n_conformers,
        ff="mmff94s",
        # ETKDG is deterministic in (molecule, seed), so pinning one seed for the
        # whole run makes every compound's geometry reproducible without needing
        # to reseed per ligand.
        etkdg_rng_seed=config.seed,
    )
    # Fallback for compounds ETKDG cannot embed: protonate but keep whatever 3D
    # coordinates the input already had. See _scrub_states().
    _TOOLS.scrub_keep_coords = Scrub(
        ph_low=config.ph,
        ph_high=None,
        skip_tautomers=config.skip_tautomers,
        skip_gen3d=True,
    )

    _TOOLS.preparator = MoleculePreparation(
        rigid_macrocycles=not config.macrocycles,
    )
    _TOOLS.config = config
    return _TOOLS


def _readable_error(exc: Exception) -> str:
    """Render an exception in terms a user can act on.

    RDKit's embedding failures arrive as a raw counter dict --
    ``{'INITIAL_COORDS': 230, 'FIRST_MINIMIZATION': 1, ...}`` -- which is the
    internal tally of which stage rejected how many attempts. That is diagnostic
    output for RDKit's authors, not something to put in a user's failure log.
    """
    text = str(exc)
    if "INITIAL_COORDS" in text or "ETK_MINIMIZATION" in text or "FINAL_CHIRAL" in text:
        return f"3D embedding failed (RDKit ETKDG could not satisfy the geometry): {text[:120]}"
    return text[:200]


def _scrub_states(tools: _Tools, mol, record: Record) -> tuple[list, str | None]:
    """Protonate and embed, falling back to the input geometry if embedding fails.

    ETKDG genuinely fails on a small fraction of strained polycyclics -- 15 in the
    first 2000 CMNPD compounds. Dropping those would be wasteful when the library
    already ships usable 3D coordinates: the compound is fine, only the attempt to
    regenerate its geometry from scratch failed.

    Returns the protomer list and a warning if the fallback was used. Falling back
    is only possible for formats that carry coordinates, so SMILES input still
    fails outright -- correctly, since there is nothing to fall back to.

    The defensive copy is essential rather than cautious: Scrubber mutates the
    molecule it is given, and strips its conformer even on the failure path. The
    fallback would otherwise be handed a molecule the first attempt had already
    emptied, and would report "no 3D conformer produced" for a compound whose
    input coordinates were perfectly good.
    """
    from rdkit import Chem

    can_fall_back = record.fmt != "smi" and tools.scrub_keep_coords is not None
    pristine = Chem.Mol(mol) if can_fall_back else None

    try:
        return list(tools.scrub(mol)), None
    except Exception as exc:  # noqa: BLE001 - any scrub failure is worth retrying
        if pristine is None:
            raise

        states = list(tools.scrub_keep_coords(pristine))
        if not states:
            raise
        return states, f"3D embedding failed ({_readable_error(exc)[:80]}); used input geometry"


def _parse(record: Record):
    """Turn a raw record into a sanitised RDKit molecule."""
    from rdkit import Chem

    if record.fmt == "smi":
        return Chem.MolFromSmiles(record.block)
    if record.fmt == "mol2":
        return Chem.MolFromMol2Block(record.block, sanitize=True, removeHs=False)
    return Chem.MolFromMolBlock(record.block, sanitize=True, removeHs=False)


def _count_torsions(pdbqt: str) -> int:
    """Count rotatable bonds in a prepared PDBQT.

    Read from the ``TORSDOF`` record, which Meeko writes at the *end* of the
    file. MGLTools' ``prepare_ligand4.py`` instead wrote an ``N active torsions``
    REMARK in the header; PDBQTs from that era are still common, so both are
    accepted, with the ``BRANCH`` count as a last resort since each branch opens
    exactly one rotatable bond.
    """
    branches = 0
    for line in pdbqt.splitlines():
        if line.startswith("TORSDOF"):
            parts = line.split()
            if len(parts) > 1 and parts[1].isdigit():
                return int(parts[1])
        elif line.startswith("BRANCH"):
            branches += 1
        elif "active torsions" in line:
            for token in line.split():
                if token.isdigit():
                    return int(token)
    return branches


def prepare_one(record: Record, config: PrepConfig, out_dir: Path) -> PreparedLigand | PrepFailure:
    """Prepare a single compound.

    Returns a failure rather than raising: one malformed molecule in a library of
    47,000 must not end the run, and the reason needs to reach the failure log
    attached to an identifier the user recognises.
    """
    from rdkit import Chem

    from drydock.core import descriptors as descmod

    try:
        mol = _parse(record)
    except Exception as exc:  # noqa: BLE001 - any RDKit failure is a prep failure
        return PrepFailure(record.ligand_id, record.compound_id, "parse", str(exc)[:200])

    if mol is None:
        return PrepFailure(
            record.ligand_id, record.compound_id, "parse", "RDKit could not parse the record"
        )

    try:
        desc = descmod.compute(mol)
    except Exception as exc:  # noqa: BLE001
        return PrepFailure(record.ligand_id, record.compound_id, "descriptors", str(exc)[:200])

    if config.max_torsions is not None and desc.rot_bonds > config.max_torsions:
        return PrepFailure(
            record.ligand_id,
            record.compound_id,
            "filter",
            f"{desc.rot_bonds} rotatable bonds exceeds max_torsions={config.max_torsions}",
        )

    tools = _get_tools(config)

    warnings: list[str] = []
    try:
        states, fallback_warning = _scrub_states(tools, mol, record)
    except Exception as exc:  # noqa: BLE001
        return PrepFailure(
            record.ligand_id, record.compound_id, "protonate", _readable_error(exc)
        )

    if fallback_warning:
        warnings.append(fallback_warning)

    if not states:
        return PrepFailure(
            record.ligand_id, record.compound_id, "protonate", "no protomer produced"
        )

    # Scrubber may return several protomers/tautomers. The first is the dominant
    # state at the requested pH, which is what should be docked; enumerating the
    # rest is a deliberate choice the user has not made here.
    prepared = states[0]

    written: list[str] = []
    torsions = 0

    out_dir.mkdir(parents=True, exist_ok=True)
    n_confs = prepared.GetNumConformers()
    if n_confs == 0:
        return PrepFailure(
            record.ligand_id, record.compound_id, "embed", "no 3D conformer produced"
        )

    try:
        from meeko import PDBQTWriterLegacy

        for conf_index in range(min(n_confs, max(1, config.n_conformers))):
            conf_id = prepared.GetConformer(conf_index).GetId()
            setups = tools.preparator.prepare(prepared, conformer_id=conf_id)
            if not setups:
                return PrepFailure(
                    record.ligand_id, record.compound_id, "pdbqt", "Meeko produced no setup"
                )

            pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
            if not ok:
                return PrepFailure(
                    record.ligand_id, record.compound_id, "pdbqt", str(err)[:200]
                )

            suffix = "" if config.n_conformers == 1 else f"_c{conf_index + 1}"
            filename = f"{record.ligand_id}{suffix}.pdbqt"
            (out_dir / filename).write_text(pdbqt, encoding="utf-8")
            written.append(filename)
            torsions = max(torsions, _count_torsions(pdbqt))
    except Exception as exc:  # noqa: BLE001
        return PrepFailure(record.ligand_id, record.compound_id, "pdbqt", str(exc)[:200])

    if torsions > TORSION_WARN_THRESHOLD:
        warnings.append(f"{torsions} active torsions; docking will be slow and imprecise")

    return PreparedLigand(
        ligand_id=record.ligand_id,
        compound_id=record.compound_id,
        pdbqt_files=tuple(written),
        torsions=torsions,
        prepared_charge=Chem.GetFormalCharge(prepared),
        descriptors=desc,
        warnings=tuple(warnings),
    )


def _worker(args: tuple[Record, PrepConfig, str]) -> PreparedLigand | PrepFailure:
    """Pool entry point.

    Embedding reproducibility is handled by the run-level ``etkdg_rng_seed`` set
    when the Scrub instance is built, not here. ETKDG is deterministic in
    ``(molecule, seed)``, so one seed for the run already gives every compound a
    reproducible geometry -- unlike docking, where the search is seeded per
    ligand because Vina explores rather than embeds.
    """
    record, config, out_dir = args
    return prepare_one(record, config, Path(out_dir))
