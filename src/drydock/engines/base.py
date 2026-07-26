"""Engine interface and shared configuration.

An engine turns one prepared ligand into a list of scored poses. Everything that
varies between AutoDock Vina, Vinardo, AutoDock4 scoring and classic AutoDock4
lives behind this interface, so the runner does not know which is in use.

Engines are constructed inside worker processes and reused across ligands. That
reuse is not an optimisation detail but the difference between a screen taking
hours and taking days -- see :class:`~drydock.engines.vina_engine.VinaEngine`.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Protocol

from drydock.core.box import Box
from drydock.core.rundir import PoseMode

# Engines Drydock can run. "ad4" is AutoDock4/AutoDock4Zn scoring driven by
# Vina's search; "autodock4" is the classic Lamarckian GA, which is 10-50x
# slower and included for small focused sets rather than library screening.
ENGINES = ("vina", "vinardo", "ad4", "autodock4")

# Engines whose scoring function reads a user-supplied parameter file. Vina's
# native function has no such file, so requesting custom parameters implies AD4.
PARAMETERISED_ENGINES = ("ad4", "autodock4")


class EngineError(RuntimeError):
    """Raised when an engine cannot dock a ligand.

    Carries the ligand identifier where known, so the runner can journal the
    failure against a name rather than an index.
    """

    def __init__(self, message: str, ligand_id: str | None = None) -> None:
        super().__init__(message)
        self.ligand_id = ligand_id


@dataclass(frozen=True, slots=True)
class DockConfig:
    """Everything needed to dock, and nothing that varies per ligand.

    Frozen and picklable so it can be handed to worker processes unchanged, and
    recorded verbatim in the run's provenance.
    """

    receptor: str
    box: Box
    engine: str = "vina"
    exhaustiveness: int = 8
    n_modes: int = 9
    energy_range: float = 3.0
    seed: int = 0

    cpu: int = 1
    """Threads per docking job.

    Deliberately 1. Vina is not reproducible across threads even with a fixed
    seed -- the search is parallel and thread scheduling perturbs which minima
    are found. Reproducibility therefore requires single-threaded jobs, with
    parallelism across ligands instead. That is also the faster arrangement for
    a library screen, so nothing is given up.
    """

    maps_dir: str | None = None
    """Directory of pre-computed AutoGrid maps, for the ad4 engine."""

    parameter_file: str | None = None
    """Custom AD4 parameter file. Implies an AD4-scoring engine."""

    def __post_init__(self) -> None:
        if self.engine not in ENGINES:
            raise ValueError(f"unknown engine {self.engine!r}; expected one of {ENGINES}")
        if self.exhaustiveness < 1:
            raise ValueError("exhaustiveness must be at least 1")
        if self.n_modes < 1:
            raise ValueError("n_modes must be at least 1")
        if self.parameter_file and self.engine not in PARAMETERISED_ENGINES:
            raise ValueError(
                f"engine {self.engine!r} has no user-editable parameter file; "
                f"custom parameters require one of {PARAMETERISED_ENGINES}"
            )

    @property
    def scoring_function(self) -> str:
        """The Vina scoring function name for this engine."""
        return "ad4" if self.engine == "ad4" else self.engine

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["box"] = self.box.to_dict()
        return data


@dataclass(frozen=True, slots=True)
class DockOutcome:
    """What an engine produces for one ligand."""

    modes: tuple[PoseMode, ...]
    poses_pdbqt: str | None = None
    """Pose coordinates, kept only for ligands worth exporting."""

    @property
    def best_affinity(self) -> float | None:
        return min((m.affinity for m in self.modes), default=None)


class Engine(Protocol):
    """What the runner requires of a docking engine."""

    name: str

    def dock(self, ligand_path: str, ligand_id: str, want_poses: bool = False) -> DockOutcome:
        """Dock one prepared ligand.

        Args:
            ligand_path: Path to the ligand PDBQT.
            ligand_id: Identifier, for error reporting.
            want_poses: Whether to return pose coordinates as well as scores.
                Skipped by default: pose text is far larger than the scores and
                only a small fraction of a screen is ever looked at.

        Raises:
            EngineError: If the ligand cannot be docked.
        """
        ...


def parse_vina_poses(pdbqt: str) -> tuple[PoseMode, ...]:
    """Extract scored modes from Vina's pose output.

    The RMSD columns come from here rather than from ``Vina.energies()``, which
    reports only energy terms. Vina writes one ``REMARK VINA RESULT`` line per
    model, carrying affinity and the two RMSDs relative to the best mode -- the
    same three numbers AutoDock Vina prints to the terminal and that PaDEL-ADV
    recorded.
    """
    modes: list[PoseMode] = []
    for line in pdbqt.splitlines():
        if "VINA RESULT" not in line:
            continue
        parts = line.split(":", 1)[-1].split()
        if len(parts) < 3:
            continue
        try:
            modes.append(
                PoseMode(
                    mode=len(modes) + 1,
                    affinity=float(parts[0]),
                    rmsd_lb=float(parts[1]),
                    rmsd_ub=float(parts[2]),
                )
            )
        except ValueError:
            continue
    return tuple(modes)
