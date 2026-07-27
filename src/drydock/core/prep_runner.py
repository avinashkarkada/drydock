"""Drives ligand preparation over a whole library.

Preparation is embarrassingly parallel, compounds do not interact. So this is
a process pool over a streaming reader, with results written to disk as they
arrive rather than collected and written at the end.

The output directory is the contract between preparation and screening:

    ligands/
    |-- pdbqt/<ligand_id>.pdbqt   one per ligand (or per conformer)
    |-- manifest.csv              descriptors, one row per ligand
    |-- failures.csv              what could not be prepared, and why
    `-- prep.json                 settings and counts, for provenance

``manifest.csv`` is what makes a hit list readable: screening produces
identifiers and affinities, and this is where the chemistry to join against
lives.
"""

from __future__ import annotations

import csv
import json
import multiprocessing as mp
import os
import time
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from drydock.core.library import Record, count_records, iter_library
from drydock.core.ligprep import PrepConfig, PrepFailure, _worker

# Explicit column order so the manifest is stable across runs and diffable.
# Identifiers first, then what was produced, then the chemistry.
MANIFEST_COLUMNS: tuple[str, ...] = (
    "ligand_id",
    "compound_id",
    "n_conformers",
    "pdbqt",
    "torsions",
    "prepared_charge",
    "smiles",
    "formula",
    "mw",
    "heavy_atoms",
    "rot_bonds",
    "clogp",
    "tpsa",
    "hbd",
    "hba",
    "formal_charge",
    "rings",
    "aromatic_rings",
    "fraction_csp3",
    "warnings",
)

FAILURE_COLUMNS: tuple[str, ...] = ("ligand_id", "compound_id", "stage", "error")

# Compounds handed to a worker at a time. Preparation takes ~40 ms each, so
# per-item dispatch would spend a meaningful fraction of the run on IPC; batching
# amortises it without making progress reporting noticeably lumpy.
CHUNK_SIZE = 64

# How often to flush the manifest. Preparation is fast enough that per-row
# flushing measurably slows it, and unlike the docking journal the cost of losing
# a few rows is small: they are simply re-prepared on resume.
FLUSH_EVERY = 200


@dataclass(slots=True)
class PrepSummary:
    """Outcome of a preparation run."""

    total: int = 0
    prepared: int = 0
    failed: int = 0
    skipped: int = 0
    elapsed_s: float = 0.0
    out_dir: str = ""

    @property
    def rate_per_s(self) -> float:
        done = self.prepared + self.failed
        return done / self.elapsed_s if self.elapsed_s > 0 else 0.0

    def to_dict(self) -> dict[str, Any]:
        from dataclasses import asdict

        return asdict(self)


class PrepDir:
    """Layout of a prepared ligand directory."""

    def __init__(self, path: str | os.PathLike[str]) -> None:
        self.path = Path(path).expanduser().resolve()

    @property
    def pdbqt_dir(self) -> Path:
        return self.path / "pdbqt"

    @property
    def manifest_file(self) -> Path:
        return self.path / "manifest.csv"

    @property
    def failures_file(self) -> Path:
        return self.path / "failures.csv"

    @property
    def info_file(self) -> Path:
        return self.path / "prep.json"

    def create(self) -> PrepDir:
        self.pdbqt_dir.mkdir(parents=True, exist_ok=True)
        return self

    def prepared_ids(self) -> set[str]:
        """Ligand IDs already in the manifest, for resuming.

        A row whose PDBQT is missing is not counted as done, that combination
        means the run died between writing the file and flushing the manifest, or
        someone deleted the structures, and in both cases the ligand should be
        prepared again.
        """
        if not self.manifest_file.exists():
            return set()

        done: set[str] = set()
        with open(self.manifest_file, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                ligand_id = row.get("ligand_id")
                if not ligand_id:
                    continue
                files = (row.get("pdbqt") or "").split(";")
                if files and all((self.pdbqt_dir / f).exists() for f in files if f):
                    done.add(ligand_id)
        return done

    def read_manifest(self) -> list[dict[str, Any]]:
        if not self.manifest_file.exists():
            return []
        with open(self.manifest_file, newline="", encoding="utf-8") as fh:
            return list(csv.DictReader(fh))

    def progress(self) -> tuple[int, int]:
        """Rows written so far, as ``(prepared, failed)``.

        Counts newlines rather than parsing CSV, so a watcher can poll a
        preparation run of any size in constant time. Approximate by a row or two
        while a flush is in flight, which is the right trade for a progress
        indicator, the manifest itself remains the authority.
        """
        return (_count_data_rows(self.manifest_file), _count_data_rows(self.failures_file))

    def is_running(self) -> bool:
        """True if preparation has started but not written its summary."""
        return self.manifest_file.exists() and not self.info_file.exists()

    def read_info(self) -> dict[str, Any] | None:
        """The summary written when preparation finished, if it has."""
        if not self.info_file.exists():
            return None
        try:
            return json.loads(self.info_file.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            return None


def _count_data_rows(path: Path) -> int:
    """Count lines in a CSV, excluding the header. Zero if absent."""
    if not path.exists():
        return 0
    try:
        with open(path, "rb") as fh:
            lines = sum(1 for _ in fh)
    except OSError:
        return 0
    return max(0, lines - 1)


def _feed(
    library: Path,
    config: PrepConfig,
    pdbqt_dir: Path,
    id_field: str | None,
    fmt: str | None,
    skip: set[str],
    limit: int | None,
) -> Iterator[tuple[Record, PrepConfig, str]]:
    """Yield work items, skipping anything already prepared.

    ``limit`` bounds how far into the library to read, not how much work to emit.
    The distinction matters on resume: counting emitted items would make
    ``--limit 2000`` mean "prepare 2000 *more*", so re-running the same command
    after an interrupted run would march through the library instead of finishing
    the job and stopping. Bounding the scan makes a repeated command idempotent.
    """
    considered = 0
    out = str(pdbqt_dir)
    for record in iter_library(library, fmt=fmt, id_field=id_field):
        if limit is not None and considered >= limit:
            return
        considered += 1
        if record.ligand_id in skip:
            continue
        yield record, config, out


def run_prep(
    library: str | os.PathLike[str],
    out_dir: str | os.PathLike[str],
    config: PrepConfig | None = None,
    *,
    n_workers: int = 0,
    id_field: str | None = None,
    fmt: str | None = None,
    limit: int | None = None,
    resume: bool = True,
    progress: Callable[[PrepSummary], None] | None = None,
) -> PrepSummary:
    """Prepare an entire compound library.

    Args:
        library: Input SDF/SMILES/MOL2, optionally gzipped.
        out_dir: Directory to populate.
        config: Preparation settings.
        n_workers: Processes to use. 0 means one per CPU.
        id_field: Property holding the compound identifier.
        fmt: Library format, inferred from the filename if omitted.
        limit: Stop after this many compounds. Useful for a trial run.
        resume: Skip ligands already present in the manifest.
        progress: Called periodically with a running summary.

    Returns:
        A :class:`PrepSummary`.
    """
    library = Path(library)
    config = config or PrepConfig()
    prep = PrepDir(out_dir).create()
    n_workers = n_workers or (os.cpu_count() or 1)

    skip = prep.prepared_ids() if resume else set()
    total_records = count_records(library, fmt=fmt)

    summary = PrepSummary(
        total=total_records if limit is None else min(limit + len(skip), total_records),
        skipped=len(skip),
        out_dir=str(prep.path),
    )

    manifest_exists = prep.manifest_file.exists() and skip
    failures_exist = prep.failures_file.exists() and skip

    started = time.perf_counter()
    with (
        open(prep.manifest_file, "a" if manifest_exists else "w", newline="", encoding="utf-8")
        as manifest_fh,
        open(prep.failures_file, "a" if failures_exist else "w", newline="", encoding="utf-8")
        as failures_fh,
    ):
        manifest = csv.DictWriter(manifest_fh, fieldnames=MANIFEST_COLUMNS, extrasaction="ignore")
        failures = csv.DictWriter(failures_fh, fieldnames=FAILURE_COLUMNS, extrasaction="ignore")
        if not manifest_exists:
            manifest.writeheader()
        if not failures_exist:
            failures.writeheader()

        work = _feed(library, config, prep.pdbqt_dir, id_field, fmt, skip, limit)

        # "forkserver", not the two more obvious alternatives.
        #
        # "spawn" re-imports __main__ in every child. A user who calls run_prep()
        # at module level in a script, which is a completely reasonable thing to
        # write, then has each worker re-execute the script and start its own
        # pool. That is a fork bomb, and it presents as the machine locking up
        # rather than as an error anyone can act on.
        #
        # "fork" (the historical Linux default) avoids that, but children inherit
        # the parent's entire address space, including whatever partially
        # initialised state RDKit and Meeko hold at the moment of forking.
        #
        # "forkserver" has neither problem: children are forked from a small,
        # clean server process, and __main__ is never re-imported.
        ctx = mp.get_context("forkserver")

        if n_workers == 1:
            results = (_worker(item) for item in work)
            pool = None
        else:
            pool = ctx.Pool(processes=n_workers)
            results = pool.imap_unordered(_worker, work, chunksize=CHUNK_SIZE)

        try:
            since_flush = 0
            for result in results:
                if isinstance(result, PrepFailure):
                    failures.writerow(result.to_row())
                    summary.failed += 1
                else:
                    manifest.writerow(result.to_row())
                    summary.prepared += 1

                since_flush += 1
                if since_flush >= FLUSH_EVERY:
                    manifest_fh.flush()
                    failures_fh.flush()
                    since_flush = 0
                    summary.elapsed_s = time.perf_counter() - started
                    if progress:
                        progress(summary)
        finally:
            if pool is not None:
                pool.terminate()
                pool.join()

    summary.elapsed_s = time.perf_counter() - started
    if progress:
        progress(summary)

    _write_info(prep, library, config, summary, n_workers)
    return summary


def _write_info(
    prep: PrepDir,
    library: Path,
    config: PrepConfig,
    summary: PrepSummary,
    n_workers: int,
) -> None:
    """Record how this directory was produced, for provenance."""
    from drydock import __version__

    info = {
        "drydock_version": __version__,
        "library": str(library.resolve()),
        "library_bytes": library.stat().st_size if library.exists() else None,
        "config": config.to_dict(),
        "n_workers": n_workers,
        "summary": summary.to_dict(),
        "created_at": time.time(),
    }
    prep.info_file.write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
