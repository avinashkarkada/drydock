"""Drives a screen: prepared ligands in, journalled results out.

The shape mirrors ligand preparation: a process pool over a stream of work,
with results written to disk as they arrive. The failure model is stricter,
though, because a screen runs for hours rather than minutes and will be
interrupted at some point.

Every finished ligand is journalled and fsync'd by the parent before the run
moves on, so a killed screen loses at most the ligand in flight. See
:mod:`drydock.core.rundir`.

Parallelism
-----------

One single-threaded docking job per worker, parallel across ligands, rather than
one multi-threaded job at a time. Two reasons, and the second is the important
one:

* Vina's own threading scales poorly past a few cores for a single ligand.
* Vina is **not reproducible across threads**, even with a fixed seed: the
  parallel search means thread scheduling decides which minima are found. Only
  single-threaded jobs give the same answer twice.

So the arrangement that is faster is also the one that is reproducible.
"""

from __future__ import annotations

import csv
import multiprocessing as mp
import os
import time
from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

from drydock.core.rundir import LigandResult, RunDir, RunStatus
from drydock.engines.base import DockConfig, EngineError


class SetupFailure(RuntimeError):
    """Raised when a run cannot be set up, rather than one ligand failing.

    Distinct from EngineError so the runner can stop the whole screen instead of
    recording the same failure once per compound.
    """

# Poses are kept only for ligands that score well enough to be worth looking at.
# Exporting all of them for a 47,000-compound screen would write gigabytes that
# nobody opens, and the run directory should stay small enough to copy.
DEFAULT_POSE_CUTOFF_RANK = 500

# Ligands per worker dispatch. One, on purpose.
#
# multiprocessing returns a chunk's results only once the whole chunk is done, so
# any chunking delays the journal by (chunk size - 1) ligands. Docking takes
# seconds to minutes per ligand, against microseconds of dispatch overhead, so
# there is nothing to amortise and real cost to batching: progress reporting goes
# lumpy, an interrupted run loses the whole in-flight chunk rather than one
# ligand, and slow ligands landing in one chunk leave a worker idle at the end.
#
# Measured against a natural-product library, chunking by 4 delayed the first
# journal entry by several minutes for no gain.
CHUNK_SIZE = 1

# How often the cached status summary is rewritten. The journal is always
# current; this only paces the convenience file the GUI polls.
STATUS_EVERY = 10

# Consecutive setup failures before a run gives up.
#
# A setup failure, missing grid maps, an unreadable receptor, affects every
# ligand equally, so continuing means failing identically tens of thousands of
# times and burying the cause. Stopping after a handful reports it once, while
# still tolerating a worker that fails to start for an unrelated transient reason.
SETUP_FAILURE_LIMIT = 5


@dataclass(frozen=True, slots=True)
class LigandJob:
    """One unit of work."""

    ligand_id: str
    path: str
    want_poses: bool = False


# Worker-process state. Building an engine means loading the receptor and
# computing grid maps (~1.5 s), so it is done once per process and reused. This
# is only correct because Vina re-seeds per dock() call. See VinaEngine.
_ENGINE = None
_ENGINE_CONFIG: DockConfig | None = None


def _get_engine(config: DockConfig):
    global _ENGINE, _ENGINE_CONFIG
    if _ENGINE is None or _ENGINE_CONFIG != config:
        from drydock.engines.vina_engine import build_engine

        _ENGINE = build_engine(config)
        _ENGINE_CONFIG = config
    return _ENGINE


def _dock_one(args: tuple[LigandJob, DockConfig]) -> tuple[LigandResult, str | None]:
    """Dock one ligand in a worker. Returns the journal record and any pose text.

    Failures are returned rather than raised: one bad ligand in tens of thousands
    must not end a screen, and the reason belongs in the journal against an
    identifier the user recognises.
    """
    job, config = args
    started = time.perf_counter()

    try:
        engine = _get_engine(config)
        outcome = engine.dock(job.path, job.ligand_id, want_poses=job.want_poses)
    except EngineError as exc:
        return (
            LigandResult(
                ligand_id=job.ligand_id,
                status="failed",
                seed=config.seed,
                elapsed_s=time.perf_counter() - started,
                error=str(exc)[:300],
                error_kind="setup" if exc.setup else "ligand",
            ),
            None,
        )
    except Exception as exc:  # noqa: BLE001 - a worker must never take the run down
        return (
            LigandResult(
                ligand_id=job.ligand_id,
                status="failed",
                seed=config.seed,
                elapsed_s=time.perf_counter() - started,
                error=f"{exc.__class__.__name__}: {str(exc)[:250]}",
                error_kind="ligand",
            ),
            None,
        )

    return (
        LigandResult(
            ligand_id=job.ligand_id,
            status="ok",
            seed=config.seed,
            elapsed_s=time.perf_counter() - started,
            modes=outcome.modes,
        ),
        outcome.poses_pdbqt,
    )


def iter_ligands(pdbqt_dir: str | os.PathLike[str]) -> Iterator[tuple[str, Path]]:
    """Yield ``(ligand_id, path)`` for every prepared ligand, in stable order.

    Sorted so that a partial run covers a predictable prefix of the library,
    which makes an interrupted screen easier to reason about than one that
    stopped at an arbitrary scattering of compounds.
    """
    for path in sorted(Path(pdbqt_dir).glob("*.pdbqt")):
        yield path.stem, path


def _feed(
    ligands: Sequence[tuple[str, Path]],
    config: DockConfig,
    done: set[str],
    pose_cutoff: int,
) -> Iterator[tuple[LigandJob, DockConfig]]:
    """Yield work items, skipping anything already journalled."""
    for index, (ligand_id, path) in enumerate(ligands):
        if ligand_id in done:
            continue
        # Poses cannot be exported by score, because the score is not known until
        # after docking. Keeping them for a leading slice of the run is a
        # rough approximation. Re-docking a handful of top hits afterwards is
        # cheap, and that is what 'drydock poses' is for.
        yield LigandJob(ligand_id, str(path), want_poses=index < pose_cutoff), config


def run_screen(
    run_dir: str | os.PathLike[str],
    pdbqt_dir: str | os.PathLike[str],
    config: DockConfig,
    *,
    n_workers: int = 0,
    limit: int | None = None,
    resume: bool = True,
    pose_cutoff: int = DEFAULT_POSE_CUTOFF_RANK,
    progress: Callable[[RunStatus], None] | None = None,
) -> RunStatus:
    """Screen a prepared ligand directory against a receptor.

    Args:
        run_dir: Directory to write the run into.
        pdbqt_dir: Prepared ligand PDBQTs.
        config: Receptor, box, engine and search settings.
        n_workers: Processes to use. 0 means one per CPU.
        limit: Only consider this many ligands. For trial runs.
        resume: Skip ligands already journalled.
        pose_cutoff: Export poses for this many ligands.
        progress: Called periodically with the running status.

    Returns:
        The final :class:`~drydock.core.rundir.RunStatus`.
    """
    run = RunDir(run_dir).create()
    n_workers = n_workers or (os.cpu_count() or 1)

    # Written before any docking starts, so an interrupted run still records what
    # it was attempting, which is exactly when the question gets asked.
    from drydock.core import provenance as prov
    from drydock.core.results import find_manifest

    manifest = find_manifest(pdbqt_dir)
    run.write_provenance(
        prov.build(
            receptor=config.receptor,
            ligand_dir=pdbqt_dir,
            config=config,
            manifest=manifest,
            extra={"n_workers": n_workers, "manifest_found": manifest is not None},
        )
    )

    ligands = list(iter_ligands(pdbqt_dir))
    if limit is not None:
        ligands = ligands[:limit]
    if not ligands:
        raise ValueError(f"no prepared ligands found in {pdbqt_dir}")

    done = run.completed_ids() if resume else set()

    status = RunStatus(
        state="running",
        total=len(ligands),
        started_at=time.time(),
        engine=config.engine,
    )
    for record in run.read_journal() if resume else ():
        if record.status == "ok":
            status.completed += 1
        elif record.status == "failed":
            status.failed += 1
        else:
            status.skipped += 1
    run.write_status(status)

    ctx = mp.get_context("forkserver")
    work = _feed(ligands, config, done, pose_cutoff)

    if n_workers == 1:
        results = (_dock_one(item) for item in work)
        pool = None
    else:
        pool = ctx.Pool(processes=n_workers)
        results = pool.imap_unordered(_dock_one, work, chunksize=CHUNK_SIZE)

    try:
        since_status = 0
        consecutive_setup_failures = 0
        for result, poses in results:
            # The parent is the journal's only writer, so there is no locking and
            # no chance of interleaved records.
            run.append(result)

            if result.status == "ok":
                status.completed += 1
                consecutive_setup_failures = 0
            else:
                status.failed += 1
                if result.is_retryable:
                    consecutive_setup_failures += 1
                    if consecutive_setup_failures >= SETUP_FAILURE_LIMIT:
                        status.state = "failed"
                        status.message = result.error
                        run.write_status(status)
                        raise SetupFailure(result.error or "run could not be set up")
                else:
                    consecutive_setup_failures = 0

            if poses:
                (run.poses_dir / f"{result.ligand_id}.pdbqt").write_text(
                    poses, encoding="utf-8"
                )

            since_status += 1
            if since_status >= STATUS_EVERY:
                since_status = 0
                run.write_status(status)
                if progress:
                    progress(status)
    except SetupFailure:
        raise
    except KeyboardInterrupt:
        status.state = "cancelled"
        status.message = "interrupted by user"
        run.write_status(status)
        raise
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()

    status.state = "finished"
    status.finished_at = time.time()
    run.write_status(status)
    if progress:
        progress(status)
    return status


def write_config(run: RunDir, config: DockConfig, pdbqt_dir: str) -> None:
    """Record the run's configuration alongside its results."""
    import tomli_w

    payload = {
        "engine": config.engine,
        "receptor": config.receptor,
        "ligands": str(pdbqt_dir),
        "exhaustiveness": config.exhaustiveness,
        "n_modes": config.n_modes,
        "seed": config.seed,
        "cpu_per_job": config.cpu,
        "box": config.box.to_dict(),
    }
    if config.maps_dir:
        payload["maps_dir"] = config.maps_dir
    if config.parameter_file:
        payload["parameter_file"] = config.parameter_file

    run.config_file.write_bytes(tomli_w.dumps(payload).encode("utf-8"))


def read_manifest_index(manifest_path: str | os.PathLike[str]) -> dict[str, dict[str, str]]:
    """Load the ligand manifest keyed by ligand_id, for joining to results."""
    path = Path(manifest_path)
    if not path.exists():
        return {}
    with open(path, newline="", encoding="utf-8") as fh:
        return {row["ligand_id"]: row for row in csv.DictReader(fh) if row.get("ligand_id")}
