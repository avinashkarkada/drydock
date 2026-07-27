"""Polls a run directory and emits what changed.

This is the whole of the GUI's relationship with a screening run: it reads
files, never writes them, and holds no handle on the process doing the work. So
the run survives the window closing, and the window survives the run crashing.

Polling rather than inotify, because run directories often live on NFS on
shared systems, where filesystem notifications are unreliable or missing
entirely. A one-second `stat` plus a short read costs little enough that the
extra robustness is worth more than the syscalls saved.
"""

from __future__ import annotations

from PySide6.QtCore import QObject, QTimer, Signal

from drydock.core.rundir import LigandResult, RunDir, RunStatus

DEFAULT_POLL_MS = 1000

# Cap on records handed to the model in one tick. A run that finished thousands
# of ligands while the GUI was closed would otherwise deliver them in a single
# batch and stall the event loop; spreading them over successive ticks keeps the
# interface responsive while it catches up.
MAX_RECORDS_PER_POLL = 2000


class RunWatcher(QObject):
    """Watches one run directory, emitting changes as Qt signals."""

    statusChanged = Signal(object)  # RunStatus
    recordsAdded = Signal(list)  # list[LigandResult]
    runFinished = Signal(object)  # RunStatus
    error = Signal(str)

    def __init__(self, poll_ms: int = DEFAULT_POLL_MS, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._run: RunDir | None = None
        self._offset = 0
        self._last_state: str | None = None
        # Records read from the journal but not yet emitted, because a single
        # tick delivered more than the model should absorb at once.
        self._backlog: list[LigandResult] = []
        self._timer = QTimer(self)
        self._timer.setInterval(poll_ms)
        self._timer.timeout.connect(self.poll)

    @property
    def run(self) -> RunDir | None:
        return self._run

    def attach(self, run_dir: str) -> bool:
        """Point the watcher at a run directory and read it from the top.

        Returns False if there is nothing there, rather than raising: attaching
        to a run that has not started yet is a normal thing for a user to do.
        """
        run = RunDir(run_dir)
        if not run.exists():
            self.error.emit(f"No run directory at {run_dir}")
            return False

        self._run = run
        self._offset = 0
        self._last_state = None
        self._backlog.clear()
        self.poll()
        self._timer.start()
        return True

    def detach(self) -> None:
        self._timer.stop()
        self._run = None
        self._offset = 0
        self._last_state = None
        self._backlog.clear()

    def poll(self) -> None:
        """Read whatever changed since last time.

        Broad exception handling is intentional. A watcher that dies because a
        file was mid-rename, or a directory briefly vanished, would defeat the
        point of decoupling it from the run; surfacing the message and trying
        again next tick is the correct response to every transient case here.
        """
        run = self._run
        if run is None:
            return

        try:
            new_records, self._offset = run.tail_journal(self._offset)
        except OSError as exc:
            self.error.emit(f"Could not read journal: {exc}")
            return

        # The journal offset only moves forward, so anything read has to be kept.
        # Queue it all and release at most MAX_RECORDS_PER_POLL per tick, which
        # bounds the work done in any single trip through the event loop.
        self._backlog.extend(new_records)
        if self._backlog:
            batch = self._backlog[:MAX_RECORDS_PER_POLL]
            del self._backlog[: len(batch)]
            self.recordsAdded.emit(batch)

        try:
            status = run.read_status()
        except OSError as exc:
            self.error.emit(f"Could not read status: {exc}")
            return

        if status is None:
            # status.json is a cache and may be absent early in a run or briefly
            # mid-rewrite. Rebuilding from the journal is always correct.
            status = run.rebuild_status()

        self.statusChanged.emit(status)

        if status.state != self._last_state:
            self._last_state = status.state
            if status.state in ("finished", "failed", "cancelled"):
                self.runFinished.emit(status)
        # Polling continues after the terminal state is reported.
        # status.json is written before the last journal records are necessarily
        # visible to a reader, and a "finished" run may still be draining a
        # backlog here.

    def load_all(self) -> list[LigandResult]:
        """Read every record currently in the journal.

        Used when attaching to a run that is already large, where replaying
        through the poll loop would be needlessly slow.
        """
        if self._run is None:
            return []
        records = list(self._run.read_journal())
        self._offset = (
            self._run.journal_file.stat().st_size if self._run.journal_file.exists() else 0
        )
        # These records are being handed to the caller directly, so anything
        # queued for incremental delivery would now be a duplicate.
        self._backlog.clear()
        return records

    def current_status(self) -> RunStatus | None:
        if self._run is None:
            return None
        return self._run.read_status() or self._run.rebuild_status()
