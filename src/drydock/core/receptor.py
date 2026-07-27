"""Reading and checking a prepared receptor.

Drydock does not prepare receptors, you point it at a PDBQT you made yourself.
It does *check* the one you give it, because receptor preparation has failure
modes that produce a perfectly valid file which then silently distorts every
result in the run.

The one that motivated this module: a receptor prepared without polar hydrogens.
AutoDock represents a hydrogen-bond donor as a heavy atom with an ``HD`` hydrogen
attached. Strip the hydrogens and the donors do not become weaker, they cease
to exist. Every backbone amide, every lysine, every tyrosine hydroxyl scores
acceptor-only. Nothing errors, every affinity is wrong in the same direction, and
the file looks fine.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

from drydock.core.box import DEFAULT_PADDING, Box

# AutoDock's polar-hydrogen type. Its presence is what makes a donor a donor.
POLAR_HYDROGEN = "HD"

# Metals AutoDock4Zn and the AD4 parameter set know how to treat specially.
METALS = frozenset({"ZN", "MG", "MN", "FE", "CA", "CU", "NI", "CO", "CD", "K", "NA_ION"})

# Tetrahedral zinc pseudo-atom, added for the AutoDock4Zn force field.
ZINC_PSEUDO = "TZ"

# Canonical AutoDock atom types, keyed by their upper-case form.
#
# AutoDock types are case-sensitive and two-letter elements use element casing:
# "Zn", not "ZN". Vina rejects the wrong case outright. And reports it as a C++
# overload-resolution error that reads like a caller bug rather than a fixable
# problem with the file. Tools disagree about this often enough to be worth
# repairing rather than just reporting. AMDock's zinc_pseudo.py, for one,
# upper-cases every metal it passes through.
CANONICAL_ATOM_TYPES: dict[str, str] = {
    t.upper(): t
    for t in (
        # Elements whose AutoDock type is mixed case.
        "Mg", "Mn", "Fe", "Zn", "Ca", "Cu", "Ni", "Co", "Cd", "Hg", "Na", "Ki",
        "Cl", "Br", "Si", "Se",
        # Single-letter and all-caps types, listed so they survive normalisation.
        "H", "HD", "HS", "C", "A", "N", "NA", "NS", "OA", "OS", "F", "P", "S",
        "SA", "I", "TZ", "G", "GA", "J", "Q",
    )
}


def normalize_atom_types(pdbqt_text: str) -> tuple[str, int]:
    """Repair atom-type casing in a PDBQT.

    Returns the corrected text and how many atoms were changed. Only the type
    field (columns 78-79) is touched; coordinates, charges and element names are
    left exactly as they were.
    """
    lines = pdbqt_text.splitlines(keepends=True)
    changed = 0

    for i, line in enumerate(lines):
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 79:
            continue
        current = line[77:79].strip()
        canonical = CANONICAL_ATOM_TYPES.get(current.upper())
        if canonical and canonical != current:
            lines[i] = f"{line[:77]}{canonical:<2}{line[79:]}"
            changed += 1

    return "".join(lines), changed


@dataclass(frozen=True, slots=True)
class ReceptorAtom:
    """One atom of a receptor PDBQT."""

    serial: int
    name: str
    residue: str
    chain: str
    residue_seq: int
    x: float
    y: float
    z: float
    charge: float
    atom_type: str

    @property
    def coordinates(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)

    @property
    def is_metal(self) -> bool:
        return self.atom_type.upper() in METALS


@dataclass(frozen=True, slots=True)
class ReceptorReport:
    """What a receptor contains, and what looks wrong with it."""

    path: str
    n_atoms: int
    atom_types: dict[str, int]
    chains: tuple[str, ...]
    residue_range: tuple[int, int] | None
    n_polar_hydrogens: int
    metals: tuple[str, ...]
    has_zinc_pseudo_atoms: bool
    problems: tuple[str, ...]
    notes: tuple[str, ...]

    @property
    def ok(self) -> bool:
        return not self.problems


def read_pdbqt(path: str | Path) -> list[ReceptorAtom]:
    """Parse a receptor PDBQT.

    Read by fixed column positions rather than by splitting on whitespace: PDB
    derivatives are column-formatted, and fields run together as soon as a
    coordinate reaches four digits or a B-factor reaches three. Whitespace
    splitting works on most files and then misparses exactly the large structures
    where a mistake is hardest to notice.
    """
    atoms: list[ReceptorAtom] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                atoms.append(
                    ReceptorAtom(
                        serial=int(line[6:11]),
                        name=line[12:16].strip(),
                        residue=line[17:20].strip(),
                        chain=line[21:22].strip(),
                        residue_seq=int(line[22:26]),
                        x=float(line[30:38]),
                        y=float(line[38:46]),
                        z=float(line[46:54]),
                        charge=float(line[70:76]) if len(line) > 70 else 0.0,
                        atom_type=line[77:79].strip() if len(line) > 77 else "",
                    )
                )
            except (ValueError, IndexError):
                # A malformed atom line is worth skipping rather than aborting;
                # inspect() reports the resulting atom count so a wholesale parse
                # failure is still obvious.
                continue
    return atoms


def inspect(path: str | Path) -> ReceptorReport:
    """Check a prepared receptor and report anything suspicious.

    Separates *problems* (things that will distort results) from *notes*, which
    are just worth knowing before committing to a long run.
    """
    path = Path(path)
    atoms = read_pdbqt(path)

    if not atoms:
        return ReceptorReport(
            path=str(path),
            n_atoms=0,
            atom_types={},
            chains=(),
            residue_range=None,
            n_polar_hydrogens=0,
            metals=(),
            has_zinc_pseudo_atoms=False,
            problems=("no ATOM or HETATM records found; is this a PDBQT?",),
            notes=(),
        )

    types: dict[str, int] = {}
    for atom in atoms:
        types[atom.atom_type] = types.get(atom.atom_type, 0) + 1

    polar_h = types.get(POLAR_HYDROGEN, 0)
    metals = tuple(sorted({a.atom_type.upper() for a in atoms if a.is_metal}))
    chains = tuple(sorted({a.chain for a in atoms if a.chain}))
    residue_seqs = [a.residue_seq for a in atoms]
    has_tz = ZINC_PSEUDO in types

    problems: list[str] = []
    notes: list[str] = []

    if polar_h == 0:
        problems.append(
            "no polar hydrogens (HD atoms). AutoDock represents a hydrogen-bond "
            "donor as a heavy atom with an HD hydrogen attached, so this receptor "
            "has no donors at all. Backbone amides, Lys, Arg and Ser/Thr/Tyr "
            "hydroxyls will all score acceptor-only. Re-prepare with polar "
            "hydrogens added."
        )
    elif polar_h < len(atoms) * 0.02:
        notes.append(
            f"only {polar_h} polar hydrogens for {len(atoms)} atoms, which is fewer "
            "than a fully protonated structure would normally carry"
        )

    if not types.get("A") and not types.get("C"):
        problems.append("no carbon atoms; this does not look like a protein receptor")

    unknown = {t for t in types if t and not t.isalnum()}
    if unknown:
        problems.append(f"unrecognised atom types: {sorted(unknown)}")

    if "" in types:
        problems.append(
            f"{types['']} atoms have no AutoDock type in columns 78-79; "
            "the file may be a PDB rather than a PDBQT"
        )

    if metals:
        notes.append(f"metals present: {', '.join(metals)}")
        if "ZN" in metals and not has_tz:
            notes.append(
                "zinc present but no TZ pseudo-atoms. Plain Vina scoring treats "
                "zinc only generically; for a metalloprotein consider the ad4 "
                "engine, which needs 'drydock add-zinc-pseudo' first."
            )

    if has_tz:
        notes.append(
            f"{types[ZINC_PSEUDO]} zinc pseudo-atoms (TZ) present; this receptor "
            "is set up for AutoDock4Zn and should be used with the ad4 engine"
        )

    return ReceptorReport(
        path=str(path),
        n_atoms=len(atoms),
        atom_types=dict(sorted(types.items(), key=lambda kv: -kv[1])),
        chains=chains,
        residue_range=(min(residue_seqs), max(residue_seqs)),
        n_polar_hydrogens=polar_h,
        metals=metals,
        has_zinc_pseudo_atoms=has_tz,
        problems=tuple(problems),
        notes=tuple(notes),
    )


def select_residues(
    atoms: Iterable[ReceptorAtom],
    residues: Sequence[int],
    chain: str | None = None,
    sidechains_only: bool = False,
) -> list[ReceptorAtom]:
    """Pick out the atoms of named residues.

    Args:
        atoms: Receptor atoms.
        residues: Residue sequence numbers.
        chain: Restrict to one chain. Necessary for multi-chain receptors, where
            the same residue number appears more than once.
        sidechains_only: Exclude backbone N/CA/C/O. A pocket is lined by side
            chains, so excluding backbone tightens the box around what actually
            forms the site.
    """
    wanted = set(residues)
    backbone = {"N", "CA", "C", "O", "OXT"}

    selected = []
    for atom in atoms:
        if atom.residue_seq not in wanted:
            continue
        if chain and atom.chain != chain:
            continue
        if sidechains_only and atom.name in backbone:
            continue
        selected.append(atom)
    return selected


def box_from_residues(
    path: str | Path,
    residues: Sequence[int],
    padding: float = DEFAULT_PADDING,
    chain: str | None = None,
    sidechains_only: bool = False,
    cubic: bool = False,
) -> tuple[Box, list[ReceptorAtom]]:
    """Build a search box enclosing the named residues.

    Returns the box and the atoms it was derived from, so a caller can report
    which residues were actually found, a mistyped residue number would
    otherwise shrink the box silently.
    """
    atoms = read_pdbqt(path)
    if not atoms:
        raise ValueError(f"no atoms read from {path}")

    selected = select_residues(atoms, residues, chain=chain, sidechains_only=sidechains_only)
    if not selected:
        present = sorted({a.residue_seq for a in atoms})
        raise ValueError(
            f"none of residues {list(residues)} were found"
            + (f" in chain {chain}" if chain else "")
            + f". The receptor covers {present[0]}-{present[-1]}. "
            "Check numbering, construct residue numbers often differ from UniProt."
        )

    found = {a.residue_seq for a in selected}
    missing = sorted(set(residues) - found)
    if missing:
        raise ValueError(
            f"residues not present in the receptor: {missing}. "
            "Check numbering, construct residue numbers often differ from UniProt."
        )

    box = Box.from_atoms([a.coordinates for a in selected], padding=padding, cubic=cubic)
    return box, selected
