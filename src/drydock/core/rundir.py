"""The run directory: the single source of truth for a screening run.

Everything a screen produces lands in one directory, and everything that reads a
screen -- the GUI, ``drydock report``, a resumed run -- reads it from there.
Nothing of consequence lives in memory, which is what lets the screening process
be detached from the process watching it.

Layout::

    run/
    |-- config.toml            the run configuration
    |-- provenance.json        versions, seed, box, input checksums
    |-- journal.jsonl          append-only, one record per finished ligand
    |-- status.json            aggregate summary, rewritten periodically
    |-- poses/<id>.pdbqt       exported poses, top-N per ligand
    |-- logs/                  engine stderr for failures worth investigating
    |-- results.csv            ranked, one row per compound (derived)
    `-- results_all_modes.csv  every pose, PaDEL-ADV schema (derived)

Two properties matter and drive the design.

**Crash safety.** ``journal.jsonl`` is append-only and written by exactly one
process -- the parent, which collects finished work from its pool. Workers never
touch it, so there is no locking and no interleaving. A run killed mid-write
leaves at most one truncated final line, which :func:`read_journal` discards. A
resumed run loses at most the ligand that was in flight.

**Cheap polling.** A watcher should not re-parse a 47,000-line journal every
second. ``status.json`` carries the aggregate counts, and is written atomically
so a reader never observes a half-written file. Watchers that want individual
records tail the journal from a byte offset via :meth:`RunDir.tail_journal`.
"""

from __future__ import annotations

import json
import os
import tempfile
import time
from collections.abc import Iterator
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Literal

# Bumped when the on-disk layout changes incompatibly. Readers refuse anything
# from the future rather than silently misinterpreting it.
SCHEMA_VERSION = 1

LigandStatus = Literal["ok", "failed", "skipped"]
RunState = Literal["pending", "running", "finished", "failed", "cancelled"]


@dataclass(frozen=True, slots=True)
class PoseMode:
    """One binding mode reported by the engine.

    Field names mirror what AutoDock Vina prints, and flow through to the
    PaDEL-ADV-compatible ``results_all_modes.csv`` unchanged.
    """

    mode: int
    affinity: float
    rmsd_lb: float = 0.0
    rmsd_ub: float = 0.0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> PoseMode:
        return cls(
            mode=int(d["mode"]),
            affinity=float(d["affinity"]),
            rmsd_lb=float(d.get("rmsd_lb", 0.0)),
            rmsd_ub=float(d.get("rmsd_ub", 0.0)),
        )


@dataclass(frozen=True, slots=True)
class LigandResult:
    """The outcome of docking one ligand: one journal line.

    A failed ligand is recorded just as deliberately as a successful one. Both
    count as "done" for resume purposes, so a compound that reliably crashes the
    engine does not get retried forever on every restart.
    """

    ligand_id: str
    status: LigandStatus
    seed: int | None = None
    elapsed_s: float = 0.0
    modes: tuple[PoseMode, ...] = ()
    error: str | None = None
    timestamp: float = field(default_factory=time.time)

    @property
    def best_affinity(self) -> float | None:
        """Affinity of the top-scoring mode, or None if there are no modes."""
        if not self.modes:
            return None
        return min(m.affinity for m in self.modes)

    def to_json(self) -> str:
        payload: dict[str, Any] = {
            "ligand_id": self.ligand_id,
            "status": self.status,
            "seed": self.seed,
            "elapsed_s": round(self.elapsed_s, 4),
            "timestamp": round(self.timestamp, 3),
        }
        if self.modes:
            payload["modes"] = [m.to_dict() for m in self.modes]
        if self.error is not None:
            payload["error"] = self.error
        return json.dumps(payload, separators=(",", ":"))

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> LigandResult:
        return cls(
            ligand_id=str(d["ligand_id"]),
            status=d.get("status", "ok"),
            seed=d.get("seed"),
            elapsed_s=float(d.get("elapsed_s", 0.0)),
            modes=tuple(PoseMode.from_dict(m) for m in d.get("modes", ())),
            error=d.get("error"),
            timestamp=float(d.get("timestamp", 0.0)),
        )


@dataclass(slots=True)
class RunStatus:
    """Aggregate progress, cheap for a watcher to poll.

    Derived entirely from the journal, so it is a cache and never authoritative.
    A corrupt or stale ``status.json`` costs a rebuild, not a run.
    """

    state: RunState = "pending"
    total: int = 0
    completed: int = 0
    failed: int = 0
    skipped: int = 0
    started_at: float | None = None
    updated_at: float | None = None
    finished_at: float | None = None
    engine: str | None = None
    message: str | None = None
    schema_version: int = SCHEMA_VERSION

    @property
    def done(self) -> int:
        """Ligands the run will not attempt again, successful or not."""
        return self.completed + self.failed + self.skipped

    @property
    def remaining(self) -> int:
        return max(0, self.total - self.done)

    @property
    def fraction(self) -> float:
        return (self.done / self.total) if self.total else 0.0

    @property
    def rate_per_s(self) -> float | None:
        """Throughput in ligands/second, or None before anything finishes."""
        if not self.started_at or self.done == 0:
            return None
        elapsed = (self.updated_at or time.time()) - self.started_at
        return (self.done / elapsed) if elapsed > 0 else None

    @property
    def eta_s(self) -> float | None:
        """Projected seconds to completion, or None if not yet estimable."""
        rate = self.rate_per_s
        if not rate or self.remaining == 0:
            return None
        return self.remaining / rate

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> RunStatus:
        known = {f for f in cls.__slots__}
        return cls(**{k: v for k, v in d.items() if k in known})


def _atomic_write(path: Path, text: str) -> None:
    """Replace ``path`` with ``text``, never leaving a partial file visible.

    Written to a sibling temp file, fsync'd, then renamed. ``os.replace`` is
    atomic within a filesystem, so a concurrent reader sees either the old file
    or the new one -- which is what makes status.json safe to poll without a
    lock.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(text)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


class RunDir:
    """Handle on a run directory.

    Safe to open concurrently from a writer (the screening process) and any
    number of readers (GUI, reporting). Readers never mutate anything.
    """

    def __init__(self, path: str | os.PathLike[str]) -> None:
        self.path = Path(path).expanduser().resolve()

    # -- layout ------------------------------------------------------------

    @property
    def config_file(self) -> Path:
        return self.path / "config.toml"

    @property
    def provenance_file(self) -> Path:
        return self.path / "provenance.json"

    @property
    def journal_file(self) -> Path:
        return self.path / "journal.jsonl"

    @property
    def status_file(self) -> Path:
        return self.path / "status.json"

    @property
    def poses_dir(self) -> Path:
        return self.path / "poses"

    @property
    def logs_dir(self) -> Path:
        return self.path / "logs"

    @property
    def results_file(self) -> Path:
        return self.path / "results.csv"

    @property
    def all_modes_file(self) -> Path:
        return self.path / "results_all_modes.csv"

    def create(self) -> RunDir:
        """Create the directory skeleton. Idempotent, so resume is a no-op."""
        for d in (self.path, self.poses_dir, self.logs_dir):
            d.mkdir(parents=True, exist_ok=True)
        return self

    def exists(self) -> bool:
        return self.path.is_dir()

    # -- journal -----------------------------------------------------------

    def append(self, result: LigandResult) -> None:
        """Append one result. Called only by the run's single writer process.

        Flushed and fsync'd per record: a `kill -9` costs the ligand in flight
        and nothing already reported. That is a real cost on spinning disks, but
        for jobs measured in seconds per ligand it is noise, and it is what makes
        the resume guarantee honest.
        """
        self.journal_file.parent.mkdir(parents=True, exist_ok=True)
        with open(self.journal_file, "a", encoding="utf-8") as fh:
            fh.write(result.to_json() + "\n")
            fh.flush()
            os.fsync(fh.fileno())

    def read_journal(self) -> Iterator[LigandResult]:
        """Yield every complete record.

        A truncated final line -- the signature of a killed writer -- is skipped
        rather than raising, because the run it belongs to is by definition one
        we are trying to recover.
        """
        if not self.journal_file.exists():
            return
        with open(self.journal_file, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    yield LigandResult.from_dict(json.loads(line))
                except (json.JSONDecodeError, KeyError, TypeError, ValueError):
                    continue

    def tail_journal(self, offset: int) -> tuple[list[LigandResult], int]:
        """Read records added since byte ``offset``.

        Returns the new records and the offset to pass next time. Only whole
        lines are consumed, so a record still being written is left for the
        following call rather than parsed half-formed. This is what lets a
        watcher follow a long run in constant time per poll.
        """
        if not self.journal_file.exists():
            return [], offset
        size = self.journal_file.stat().st_size
        if size < offset:
            # Truncated or replaced underneath us; restart from the top.
            offset = 0
        if size == offset:
            return [], offset
        with open(self.journal_file, "rb") as fh:
            fh.seek(offset)
            chunk = fh.read(size - offset)
        last_nl = chunk.rfind(b"\n")
        if last_nl == -1:
            return [], offset
        complete, consumed = chunk[: last_nl + 1], last_nl + 1
        results: list[LigandResult] = []
        for raw in complete.decode("utf-8", errors="replace").splitlines():
            raw = raw.strip()
            if not raw:
                continue
            try:
                results.append(LigandResult.from_dict(json.loads(raw)))
            except (json.JSONDecodeError, KeyError, TypeError, ValueError):
                continue
        return results, offset + consumed

    def completed_ids(self) -> set[str]:
        """Ligand IDs the run should not attempt again.

        Includes failures deliberately -- see :class:`LigandResult`.
        """
        return {r.ligand_id for r in self.read_journal()}

    # -- status ------------------------------------------------------------

    def write_status(self, status: RunStatus) -> None:
        status.updated_at = time.time()
        _atomic_write(self.status_file, json.dumps(status.to_dict(), indent=2) + "\n")

    def read_status(self) -> RunStatus | None:
        """Read the cached status, or None if absent or unreadable.

        Returning None on a malformed file is deliberate: status.json is a cache,
        and a watcher should fall back to rebuilding from the journal rather than
        surfacing an error for something that is not authoritative.
        """
        if not self.status_file.exists():
            return None
        try:
            return RunStatus.from_dict(json.loads(self.status_file.read_text(encoding="utf-8")))
        except (json.JSONDecodeError, OSError, TypeError, ValueError):
            return None

    def rebuild_status(self, total: int | None = None) -> RunStatus:
        """Recompute status from the journal, ignoring status.json entirely."""
        cached = self.read_status()
        status = RunStatus(
            state=cached.state if cached else "pending",
            total=total if total is not None else (cached.total if cached else 0),
            started_at=cached.started_at if cached else None,
            engine=cached.engine if cached else None,
        )
        for r in self.read_journal():
            if r.status == "ok":
                status.completed += 1
            elif r.status == "failed":
                status.failed += 1
            else:
                status.skipped += 1
        status.updated_at = time.time()
        return status

    # -- provenance --------------------------------------------------------

    def write_provenance(self, provenance: dict[str, Any]) -> None:
        payload = {"schema_version": SCHEMA_VERSION, **provenance}
        _atomic_write(self.provenance_file, json.dumps(payload, indent=2) + "\n")

    def read_provenance(self) -> dict[str, Any] | None:
        if not self.provenance_file.exists():
            return None
        try:
            return json.loads(self.provenance_file.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            return None

    def __repr__(self) -> str:
        return f"RunDir({str(self.path)!r})"
