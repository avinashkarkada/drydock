"""Receptor preparation: structure in, docking-ready PDBQT out.

Originally out of scope. The reason was practical rather than principled --
receptor preparation meant MGLTools, which is Python 2, unmaintained, and a
dependency no reproducible environment should carry. Asking users to prepare
receptors elsewhere was better than pinning a dead toolchain.

Meeko's ``mk_prepare_receptor`` removes that objection: it is already in the
pinned environment, maintained by the same group as the docking engine, and adds
polar hydrogens correctly. So preparation now happens here.

What this adds over calling Meeko directly
------------------------------------------

* **mmCIF input.** Meeko's PDB reader does not take mmCIF without ProDy, so
  structures are converted with gemmi first. Most people download mmCIF now.
* **Checking the output.** Preparation can succeed and still produce a receptor
  that ruins a screen. The result is inspected and reported, so a file with no
  hydrogen-bond donors is caught here rather than after a two-day run.
* **Explaining failures.** Meeko refuses structures with residues missing atoms,
  which is correct but arrives as a traceback. That is translated into what it
  means and which flag addresses it.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from drydock.core.receptor import ReceptorReport, inspect

# Structures Meeko reads directly. Anything else is converted to PDB first.
NATIVE_FORMATS = {".pdb", ".pdbqt", ".pqr"}
CIF_FORMATS = {".cif", ".mmcif", ".bcif"}


class ReceptorPrepError(RuntimeError):
    """Raised when a receptor cannot be prepared."""


@dataclass(frozen=True, slots=True)
class PrepResult:
    """A prepared receptor and what checking it found."""

    receptor_pdbqt: Path
    report: ReceptorReport
    converted_from: Path | None = None
    deleted_residues: tuple[str, ...] = ()

    @property
    def ok(self) -> bool:
        return self.report.ok


def to_pdb(structure: str | Path, output: str | Path | None = None) -> Path:
    """Convert a structure to PDB, for formats Meeko cannot read directly.

    Uses gemmi rather than Open Babel: gemmi is an mmCIF library first, so it
    preserves chain and entity information that a general format converter
    routinely mangles.
    """
    structure = Path(structure)
    if structure.suffix.lower() in NATIVE_FORMATS:
        return structure

    if structure.suffix.lower() not in CIF_FORMATS:
        raise ReceptorPrepError(
            f"unsupported structure format {structure.suffix!r}; "
            f"expected one of {sorted(NATIVE_FORMATS | CIF_FORMATS)}"
        )

    try:
        import gemmi
    except ImportError as exc:  # pragma: no cover - gemmi is a pinned dependency
        raise ReceptorPrepError("gemmi is needed to read mmCIF files") from exc

    output = Path(output) if output else structure.with_suffix(".pdb")
    try:
        parsed = gemmi.read_structure(str(structure))
        parsed.setup_entities()
        parsed.write_pdb(str(output))
    except Exception as exc:  # noqa: BLE001 - any gemmi failure is a prep failure
        raise ReceptorPrepError(f"could not convert {structure.name}: {exc}") from exc

    return output


def prepare_receptor(
    structure: str | Path,
    output_basename: str | Path,
    *,
    allow_bad_residues: bool = True,
    delete_residues: str | None = None,
    charge_model: str = "gasteiger",
    default_altloc: str | None = None,
    keep_pdb: bool = False,
) -> PrepResult:
    """Prepare a receptor for docking.

    Args:
        structure: PDB or mmCIF input.
        output_basename: Output path without extension; ``.pdbqt`` is appended.
        allow_bad_residues: Delete residues with missing atoms rather than
            failing. On by default because crystal structures routinely have
            disordered side chains, and refusing the whole structure over one
            unresolved lysine helps nobody. Deletions are reported.
        delete_residues: Residues to remove, e.g. ``"A:350,B:15,16"``. Use for
            waters, cryoprotectants and other crystallisation artefacts.
        charge_model: ``gasteiger`` (default), ``espaloma`` or ``zero``.
        default_altloc: Which alternate location to take where a residue has
            several. Meeko fails on ambiguity rather than guessing.
        keep_pdb: Keep the intermediate PDB when converting from mmCIF.

    Returns:
        A :class:`PrepResult`, including the check of what was produced.

    Raises:
        ReceptorPrepError: If preparation fails.
    """
    structure = Path(structure)
    if not structure.exists():
        raise ReceptorPrepError(f"structure not found: {structure}")

    output_basename = Path(output_basename)
    output_basename.parent.mkdir(parents=True, exist_ok=True)

    converted_from = None
    pdb_input = structure
    if structure.suffix.lower() in CIF_FORMATS:
        pdb_input = to_pdb(structure, output_basename.with_suffix(".converted.pdb"))
        converted_from = structure

    executable = shutil.which("mk_prepare_receptor.py") or shutil.which("mk_prepare_receptor")
    if not executable:
        raise ReceptorPrepError(
            "mk_prepare_receptor.py not found. It ships with Meeko, which is part "
            "of Drydock's pinned environment -- run through 'pixi run'."
        )

    command = [
        sys.executable if executable.endswith(".py") else executable,
        *([executable] if executable.endswith(".py") else []),
        "--read_pdb", str(pdb_input),
        "-o", str(output_basename),
        "-p",
        "--charge_model", charge_model,
    ]
    if allow_bad_residues:
        command.append("-a")
    if delete_residues:
        command += ["-d", delete_residues]
    if default_altloc:
        command += ["--default_altloc", default_altloc]

    result = subprocess.run(command, capture_output=True, text=True)

    if not keep_pdb and converted_from and pdb_input.exists():
        pdb_input.unlink(missing_ok=True)

    produced = output_basename.with_suffix(".pdbqt")
    if result.returncode != 0 or not produced.exists():
        raise ReceptorPrepError(_explain_failure(result.stdout + result.stderr))

    return PrepResult(
        receptor_pdbqt=produced,
        report=inspect(produced),
        converted_from=converted_from,
        deleted_residues=_deleted_residues(result.stdout + result.stderr),
    )


def _deleted_residues(output: str) -> tuple[str, ...]:
    """Pull the residues Meeko removed out of its output.

    Worth surfacing: a silently deleted active-site residue changes the pocket,
    and the user should hear about it from the tool rather than discover it in a
    pose.
    """
    deleted = []
    for line in output.splitlines():
        lowered = line.lower()
        if "deleted" in lowered or "removing residue" in lowered:
            deleted.append(line.strip()[:120])
    return tuple(deleted)


def _explain_failure(output: str) -> str:
    """Translate Meeko's failure into something actionable."""
    text = output.strip()

    if "missing atom" in text.lower() or "incomplete" in text.lower():
        return (
            "the structure has residues with missing atoms, which Meeko refuses "
            "rather than guess at. Pass --allow-bad-residues to delete them "
            "(they are reported), or repair the structure first.\n\n"
            + text[-600:]
        )
    if "altloc" in text.lower():
        return (
            "the structure has alternate locations Meeko will not choose between. "
            "Pass --altloc A (or B) to pick one.\n\n" + text[-600:]
        )
    if "template" in text.lower():
        return (
            "a residue does not match any known template -- commonly a modified "
            "residue, a ligand left in the file, or a non-standard cofactor. "
            "Delete it with --delete-residues, or supply a template.\n\n"
            + text[-600:]
        )
    return f"mk_prepare_receptor failed:\n{text[-800:]}"
