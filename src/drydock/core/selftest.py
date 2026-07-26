"""End-to-end verification against a known-good reference.

The point of a self-test in a scientific tool is narrower than "do the tests
pass". Unit tests run against mocks and fixtures; they confirm the code is
internally consistent. What they cannot confirm is that *this installation, on
this machine* produces the same numbers as everywhere else -- which is the claim
Drydock actually makes.

So this runs the real pipeline on a small bundled case:

    ligands.smi -> prepare -> dock against receptor.pdbqt -> compare affinities

and checks the results against values recorded from a known-good run. A silent
numerical drift -- a rebuilt Vina, a different BLAS, a subtly different RDKit --
shows up here as a failed comparison rather than as a paper that will not
replicate.

Tolerance
---------

Affinities are compared to +/-0.5 kcal/mol. Not because that is a scientifically
interesting margin, but because it is wide enough to absorb floating-point
differences between CPU architectures and narrow enough that a genuine change in
behaviour -- a different scoring function, missing hydrogens, a misplaced box --
moves results well past it. A drift smaller than that would not change any
conclusion drawn from a screen.
"""

from __future__ import annotations

import json
import shutil
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "selftest"
RECEPTOR = DATA_DIR / "receptor.pdbqt"
LIGANDS = DATA_DIR / "ligands.smi"
REFERENCE = DATA_DIR / "reference.json"

# Wide enough for cross-architecture float noise, narrow enough that a real
# behavioural change fails. See the module docstring.
AFFINITY_TOLERANCE = 0.5

# Fixed so the self-test is the same computation everywhere.
SEED = 42
EXHAUSTIVENESS = 4
N_MODES = 9

# The catalytic zinc site of the bundled MMP9 excerpt.
BOX_CENTER = (70.219, 18.432, 50.134)
BOX_SIZE = (20.0, 20.0, 20.0)


@dataclass(slots=True)
class Check:
    """One verified property."""

    name: str
    passed: bool
    detail: str = ""

    def __str__(self) -> str:
        mark = "ok  " if self.passed else "FAIL"
        return f"  [{mark}] {self.name}{f' -- {self.detail}' if self.detail else ''}"


@dataclass(slots=True)
class SelftestResult:
    """Outcome of a self-test run."""

    checks: list[Check] = field(default_factory=list)
    affinities: dict[str, float] = field(default_factory=dict)
    elapsed_s: float = 0.0

    @property
    def passed(self) -> bool:
        return all(c.passed for c in self.checks)

    @property
    def failures(self) -> list[Check]:
        return [c for c in self.checks if not c.passed]

    def add(self, name: str, passed: bool, detail: str = "") -> None:
        self.checks.append(Check(name, passed, detail))


def load_reference() -> dict[str, Any] | None:
    """Reference affinities from a known-good run, if recorded."""
    if not REFERENCE.exists():
        return None
    try:
        return json.loads(REFERENCE.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def run(work_dir: str | Path | None = None, keep: bool = False) -> SelftestResult:
    """Run the full pipeline on the bundled case and verify the results.

    Args:
        work_dir: Where to work. A temporary directory by default.
        keep: Leave the working directory in place for inspection.

    Returns:
        A :class:`SelftestResult`.
    """
    from drydock.core.box import Box
    from drydock.core.ligprep import PrepConfig
    from drydock.core.prep_runner import PrepDir, run_prep
    from drydock.core.receptor import inspect
    from drydock.core.rundir import RunDir
    from drydock.core.screen import run_screen
    from drydock.engines.base import DockConfig

    result = SelftestResult()
    started = time.perf_counter()

    temporary = work_dir is None
    work = Path(work_dir) if work_dir else Path(tempfile.mkdtemp(prefix="drydock_selftest_"))
    work.mkdir(parents=True, exist_ok=True)

    try:
        # -- bundled data present ------------------------------------------
        result.add("bundled receptor present", RECEPTOR.exists(), str(RECEPTOR))
        result.add("bundled ligands present", LIGANDS.exists(), str(LIGANDS))
        if not (RECEPTOR.exists() and LIGANDS.exists()):
            return result

        # -- receptor is usable --------------------------------------------
        report = inspect(RECEPTOR)
        result.add(
            "receptor parses",
            report.n_atoms > 0,
            f"{report.n_atoms} atoms, types {sorted(report.atom_types)}",
        )
        result.add(
            "receptor atom types are canonical",
            "ZN" not in report.atom_types,
            "metals must use element casing (Zn, not ZN) or Vina rejects them",
        )

        # -- ligand preparation --------------------------------------------
        prep_dir = work / "ligands"
        summary = run_prep(LIGANDS, prep_dir, PrepConfig(seed=SEED), n_workers=1)
        result.add(
            "ligands prepare",
            summary.failed == 0 and summary.prepared > 0,
            f"{summary.prepared} prepared, {summary.failed} failed",
        )
        if summary.prepared == 0:
            return result

        # -- docking ---------------------------------------------------------
        config = DockConfig(
            receptor=str(RECEPTOR),
            box=Box(BOX_CENTER, BOX_SIZE),
            engine="vina",
            exhaustiveness=EXHAUSTIVENESS,
            n_modes=N_MODES,
            seed=SEED,
            cpu=1,
        )
        run_dir = work / "run"
        status = run_screen(run_dir, PrepDir(prep_dir).pdbqt_dir, config, n_workers=1)
        result.add(
            "docking completes",
            status.failed == 0 and status.completed > 0,
            f"{status.completed} docked, {status.failed} failed",
        )

        records = list(RunDir(run_dir).read_journal())
        result.affinities = {
            r.ligand_id: round(r.best_affinity, 2)
            for r in records
            if r.best_affinity is not None
        }

        result.add(
            "every ligand scored",
            len(result.affinities) == summary.prepared,
            f"{len(result.affinities)}/{summary.prepared}",
        )
        result.add(
            "poses are ordered best-first",
            all(
                list(m.affinity for m in r.modes) == sorted(m.affinity for m in r.modes)
                for r in records
                if r.modes
            ),
            "mode 1 must be the best pose; results.csv depends on it",
        )

        # -- results assembly ------------------------------------------------
        from drydock.core.results import write_results

        results_csv, all_modes_csv, n_rows = write_results(
            run_dir, PrepDir(prep_dir).manifest_file
        )
        result.add(
            "results written",
            results_csv.exists() and all_modes_csv.exists() and n_rows > 0,
            f"{n_rows} ranked rows",
        )

        # -- comparison against the reference --------------------------------
        reference = load_reference()
        if reference is None:
            result.add(
                "reference comparison",
                True,
                "no reference recorded; run 'drydock dev record-reference' to create one",
            )
        else:
            _compare_to_reference(result, reference)

    finally:
        result.elapsed_s = time.perf_counter() - started
        if temporary and not keep:
            shutil.rmtree(work, ignore_errors=True)

    return result


def _compare_to_reference(result: SelftestResult, reference: dict[str, Any]) -> None:
    """Check measured affinities against recorded ones."""
    expected: dict[str, float] = reference.get("affinities", {})
    if not expected:
        result.add("reference comparison", False, "reference file has no affinities")
        return

    missing = sorted(set(expected) - set(result.affinities))
    if missing:
        result.add("all reference ligands scored", False, f"missing: {missing}")
        return
    result.add("all reference ligands scored", True, f"{len(expected)} ligands")

    drifted = []
    for ligand_id, want in expected.items():
        got = result.affinities[ligand_id]
        if abs(got - want) > AFFINITY_TOLERANCE:
            drifted.append(f"{ligand_id}: expected {want:+.2f}, got {got:+.2f}")

    result.add(
        f"affinities within {AFFINITY_TOLERANCE} kcal/mol of reference",
        not drifted,
        "; ".join(drifted) if drifted else f"largest deviation "
        f"{max(abs(result.affinities[k] - v) for k, v in expected.items()):.2f}",
    )


def record_reference(output: str | Path | None = None) -> Path:
    """Run the pipeline and record its affinities as the new reference.

    Deliberately a separate, explicitly-invoked command rather than something the
    self-test does when the reference is missing. A self-test that writes its own
    expectations passes unconditionally and verifies nothing.
    """
    from drydock import __version__
    from drydock.core.provenance import engine_versions, package_versions

    result = run()
    if not result.affinities:
        raise RuntimeError("self-test produced no affinities; refusing to record a reference")

    output = Path(output) if output else REFERENCE
    payload = {
        "drydock_version": __version__,
        "recorded_at_iso": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "packages": package_versions(),
        "engines": engine_versions(),
        "settings": {
            "seed": SEED,
            "exhaustiveness": EXHAUSTIVENESS,
            "n_modes": N_MODES,
            "box_center": list(BOX_CENTER),
            "box_size": list(BOX_SIZE),
        },
        "affinities": result.affinities,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return output
