"""AutoDock4Zn support: zinc pseudo-atoms and grid parameter files.

Docking to a zinc metalloprotein with a generic scoring function is a known
weakness. Vina types metals as unremarkable hydrogen-bond donors, so it has no
notion of the tetrahedral coordination geometry that actually governs how an
inhibitor engages a catalytic zinc.

AutoDock4Zn handles this by adding pseudo-atoms (type ``TZ``) at the vacant
tetrahedral positions around each zinc, together with a set of pairwise
potentials. A ligand atom reaching a TZ site is rewarded for completing the
coordination sphere in the right geometry, not just for being nearby.

Two things are needed, both provided here:

1. **TZ atoms in the receptor**, added by :func:`add_zinc_pseudo_atoms`. The
   AD4Zn path cannot run without them.
2. **A grid parameter file** carrying the AD4Zn potentials, which AutoGrid reads
   when computing maps.

Provenance
----------

``zinc_pseudo.py`` and the TZ parameter line are vendored from AMDock (GPL-3),
which carries them from the AutoDock4Zn authors. AMDock also ships
``prepare_gpf4zn.py``, but that imports MolKit and AutoDockTools from MGLTools,
which is Python 2 and no longer practical to install. The GPF is written here
instead. Not much is lost: the AD4Zn-specific part of a GPF is the six
``nbp_r_eps`` lines below.
"""

from __future__ import annotations

import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from drydock.core.box import Box

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

ZINC_PSEUDO_SCRIPT = DATA_DIR / "zinc_pseudo.py"
AD4ZN_PARAMETERS = DATA_DIR / "AD4Zn.dat"

# Grid spacing in Angstroms. This is AutoGrid's default and what the AD4
# scoring function was parameterised against. Changing it invalidates the
# potentials.
GRID_SPACING = 0.375

# The AutoDock4Zn potentials, appended to an otherwise standard GPF.
#
# This is the entire AD4Zn-specific content of the grid parameter file. The first
# line is the reward for a ligand acceptor (NA) reaching a tetrahedral zinc
# pseudo-atom; the rest retune the direct zinc-to-ligand interactions so that the
# pseudo-atoms, rather than raw electrostatics, carry the coordination geometry.
#
# Values reproduced from AMDock's prepare_gpf4zn.py (GPL-3), which carries them
# from Santos-Martins et al., J Chem Inf Model 54:2371 (2014).
AD4ZN_POTENTIALS: tuple[str, ...] = (
    "nbp_r_eps 0.25 3.8581 12 6 NA TZ",
    "nbp_r_eps 2.1  0.6391 12 6 OA Zn",
    "nbp_r_eps 2.25 1.2617 12 6 SA Zn",
    "nbp_r_eps 1.0  0.0    12 6 HD Zn",
    "nbp_r_eps 2.0  0.001  12 6 NA Zn",
    "nbp_r_eps 2.0  0.0493 12 6  N Zn",
)

# Ligand atom types to compute maps for. AutoGrid needs a map per type, and a
# screen must cover every type any ligand might present, so this is the union of
# what Meeko emits rather than what one ligand happens to use.
LIGAND_ATOM_TYPES: tuple[str, ...] = (
    "HD", "C", "A", "N", "NA", "OA", "F", "P", "SA", "S", "Cl", "Br", "I", "Si",
)


class ZincError(RuntimeError):
    """Raised when zinc pseudo-atom placement or map generation fails."""


# Vina loads AD4 maps by this prefix, so the files are always named receptor.*
# regardless of what the input structure was called.
MAP_PREFIX = "receptor"


def maps_status(maps_dir: str | Path | None) -> tuple[bool, str]:
    """Check whether a directory holds usable AutoGrid maps.

    Without this check, Vina reports a missing map set as ``Cannot find affinity
    maps with <path>`` once per ligand, so a screen pointed at the wrong
    directory fails tens of thousands of times without ever saying what is
    actually wrong.

    The usual mistake is pointing at the directory holding the receptor. Maps do
    not live there until ``drydock maps`` has been run.

    Returns:
        ``(ok, message)``. The message says what to do when not ok.
    """
    if not maps_dir:
        return False, "no maps directory given; the ad4 engine needs one"

    path = Path(maps_dir)
    if not path.is_dir():
        return False, f"{path} is not a directory"

    field = path / f"{MAP_PREFIX}.maps.fld"
    maps = sorted(path.glob(f"{MAP_PREFIX}.*.map"))

    if field.exists() and maps:
        return True, f"{len(maps)} maps"

    receptors = sorted(path.glob("*.pdbqt"))
    if receptors and not maps:
        names = ", ".join(p.name for p in receptors[:3])
        return False, (
            f"{path} contains receptor files ({names}) but no AutoGrid maps. "
            "Maps are computed separately from the receptor, run 'drydock maps' "
            "(or use Generate maps in the GUI) and point this at the directory it "
            "writes."
        )

    if maps and not field.exists():
        return False, (
            f"{path} has {len(maps)} map files but no {field.name}. The map set is "
            "incomplete; regenerate it."
        )

    return False, (
        f"{path} contains no AutoGrid maps. Run 'drydock maps' first, or choose "
        "the vina or vinardo engine, which compute their own maps."
    )


@dataclass(frozen=True, slots=True)
class ZincResult:
    """Outcome of adding zinc pseudo-atoms to a receptor."""

    receptor_tz: Path
    n_pseudo_atoms: int
    n_zinc: int


def add_zinc_pseudo_atoms(receptor: str | Path, output: str | Path | None = None) -> ZincResult:
    """Place tetrahedral TZ pseudo-atoms around each zinc in a receptor.

    Runs the vendored ``zinc_pseudo.py``, which inspects each zinc's existing
    coordination and places pseudo-atoms at the vacant tetrahedral positions,
    2.0 A out. It is run as a subprocess rather than imported: it is third-party
    code with a script-shaped interface.

    Args:
        receptor: Prepared receptor PDBQT containing at least one zinc.
        output: Where to write. Defaults to ``<receptor stem>_TZ.pdbqt``.

    Returns:
        A :class:`ZincResult`.

    Raises:
        ZincError: If the receptor has no zinc, or placement fails.
    """
    receptor = Path(receptor)
    if not receptor.exists():
        raise ZincError(f"receptor not found: {receptor}")

    if not ZINC_PSEUDO_SCRIPT.exists():
        raise ZincError(f"vendored zinc_pseudo.py is missing from {DATA_DIR}")

    text = receptor.read_text(encoding="utf-8", errors="replace")
    n_zinc = sum(
        1
        for line in text.splitlines()
        if line.startswith(("ATOM", "HETATM")) and line[77:79].strip().upper() == "ZN"
    )
    if n_zinc == 0:
        raise ZincError(
            f"{receptor} contains no zinc atoms, so AutoDock4Zn has nothing to do. "
            "Use the vina or vinardo engine instead."
        )

    output = Path(output) if output else receptor.with_name(f"{receptor.stem}_TZ.pdbqt")

    result = subprocess.run(
        [sys.executable, str(ZINC_PSEUDO_SCRIPT), "-r", str(receptor), "-o", str(output)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0 or not output.exists():
        raise ZincError(
            f"zinc_pseudo.py failed: {(result.stderr or result.stdout).strip()[:300]}"
        )

    # zinc_pseudo.py upper-cases the metal types it passes through, turning "Zn"
    # into "ZN", which AutoGrid tolerates but Vina rejects outright. Repair the
    # casing so the output is usable with every engine, not just the ad4 path it
    # was produced for.
    from drydock.core.receptor import normalize_atom_types

    produced = output.read_text(encoding="utf-8", errors="replace")
    repaired, n_fixed = normalize_atom_types(produced)
    if n_fixed:
        output.write_text(repaired, encoding="utf-8")
        produced = repaired

    n_tz = sum(
        1
        for line in produced.splitlines()
        if line.startswith(("ATOM", "HETATM")) and line[77:79].strip() == "TZ"
    )
    if n_tz == 0:
        raise ZincError(
            "zinc_pseudo.py produced no TZ atoms. This happens when every zinc is "
            "already fully coordinated by the protein, leaving no vacant site for a "
            "ligand, in which case AutoDock4Zn has nothing to contribute."
        )

    return ZincResult(receptor_tz=output, n_pseudo_atoms=n_tz, n_zinc=n_zinc)


def write_gpf(
    receptor_tz: str | Path,
    box: Box,
    output: str | Path,
    ligand_types: tuple[str, ...] = LIGAND_ATOM_TYPES,
    parameter_file: str | Path | None = None,
    spacing: float = GRID_SPACING,
) -> Path:
    """Write an AutoGrid parameter file with the AutoDock4Zn potentials.

    Replaces AMDock's ``prepare_gpf4zn.py``, which cannot run without MGLTools.

    Args:
        receptor_tz: Receptor PDBQT including TZ pseudo-atoms.
        box: Search box. Determines grid centre and extent.
        output: Where to write the GPF.
        ligand_types: Atom types to compute maps for. Must cover every type any
            ligand in the screen presents, a missing map is a hard failure at
            docking time, not a silently skipped atom.
        parameter_file: AD4 parameter file. Defaults to the bundled AD4Zn.dat.
        spacing: Grid spacing in Angstroms.

    Returns:
        The path written.
    """
    receptor_tz = Path(receptor_tz)
    output = Path(output)
    parameters = Path(parameter_file) if parameter_file else AD4ZN_PARAMETERS

    # AutoGrid wants an odd number of points per axis so the grid is centred on a
    # point rather than between two.
    npts = []
    for size in box.size:
        n = int(round(size / spacing))
        npts.append(n if n % 2 == 0 else n + 1)

    receptor_types = _receptor_atom_types(receptor_tz)

    lines = [
        f"parameter_file {parameters}",
        f"npts {npts[0]} {npts[1]} {npts[2]}",
        "gridfld receptor.maps.fld",
        f"spacing {spacing}",
        f"receptor_types {' '.join(receptor_types)}",
        f"ligand_types {' '.join(ligand_types)}",
        f"receptor {receptor_tz.name}",
        f"gridcenter {box.center[0]:.3f} {box.center[1]:.3f} {box.center[2]:.3f}",
        "smooth 0.5",
    ]
    lines += [f"map receptor.{t}.map" for t in ligand_types]
    lines += [
        "elecmap receptor.e.map",
        "dsolvmap receptor.d.map",
        "dielectric -0.1465",
    ]
    lines += list(AD4ZN_POTENTIALS)

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output


def _receptor_atom_types(receptor: Path) -> list[str]:
    """Collect the distinct AutoDock types present in a receptor, in file order."""
    seen: list[str] = []
    for line in receptor.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom_type = line[77:79].strip()
        if atom_type and atom_type not in seen:
            seen.append(atom_type)
    return seen


def run_autogrid(gpf: str | Path, work_dir: str | Path | None = None) -> Path:
    """Run AutoGrid to compute the maps described by a GPF.

    AutoGrid resolves the paths inside a GPF relative to the working directory,
    so it is run from the directory containing the file.

    Returns:
        The directory containing the maps.
    """
    gpf = Path(gpf).resolve()
    work_dir = Path(work_dir).resolve() if work_dir else gpf.parent
    work_dir.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(
        ["autogrid4", "-p", gpf.name, "-l", f"{gpf.stem}.glg"],
        cwd=work_dir,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        log = work_dir / f"{gpf.stem}.glg"
        detail = log.read_text(errors="replace")[-800:] if log.exists() else result.stderr
        raise ZincError(f"autogrid4 failed:\n{detail.strip()[:800]}")

    if not (work_dir / "receptor.maps.fld").exists():
        raise ZincError(f"autogrid4 reported success but wrote no maps in {work_dir}")

    return work_dir
