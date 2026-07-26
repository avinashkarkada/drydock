"""Drydock's main window.

A monitor, not a controller. It reads a run directory and displays it; the
screening process is launched detached and is not owned by this window. Closing
it does not stop a run, and a run crashing does not take the window with it.

The interface follows the order the work is done:

* **1. Ligands** -- convert, protonate, optimise and prepare a library.
* **2. Screen** -- receptor, box and engine; starts a detached screen.
* **3. Results** -- a virtualised table of every ligand scored so far.

with a run bar and a bounded activity log framing them.

Starting a screen from Setup spawns ``drydock screen`` in its own session and
then simply attaches to the run directory it creates. The window holds no handle
on the work, which is what lets it be closed and reopened freely.
"""

from __future__ import annotations

import sys
from collections import deque
from pathlib import Path

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QAction, QFont
from PySide6.QtWidgets import (
    QApplication,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMainWindow,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTableView,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from drydock import __version__
from drydock.core.rundir import LigandResult, RunStatus
from drydock.gui.ligand_panel import LigandPanel
from drydock.gui.model import ResultsTableModel
from drydock.gui.setup_panel import SetupPanel
from drydock.gui.watcher import RunWatcher

# The activity pane is a rolling window, not a transcript. A long screen emits
# tens of thousands of events, and retaining them would grow the process without
# bound for information nobody scrolls back to. The engine's own output is kept
# on disk under the run's logs/ directory.
ACTIVITY_LINES = 500

# Index into ResultsTableModel.COLUMNS of the affinity column, which is the one
# worth sorting by on open: the point of a screen is the best binders.
AFFINITY_COLUMN = 2


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle(f"Drydock {__version__}")
        self.resize(1100, 720)

        self._watcher = RunWatcher(parent=self)
        self._model = ResultsTableModel(self)
        self._activity: deque[str] = deque(maxlen=ACTIVITY_LINES)

        self._build_ui()
        self._connect()

    # -- construction ------------------------------------------------------

    def _build_ui(self) -> None:
        central = QWidget()
        layout = QVBoxLayout(central)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)

        layout.addLayout(self._build_run_bar())
        layout.addWidget(self._build_progress())

        splitter = QSplitter(Qt.Orientation.Vertical)
        splitter.addWidget(self._build_tabs())
        splitter.addWidget(self._build_activity())
        splitter.setStretchFactor(0, 4)
        splitter.setStretchFactor(1, 1)
        layout.addWidget(splitter, stretch=1)

        self.setCentralWidget(central)
        self.statusBar().showMessage("No run attached")
        self._build_menu()

    def _build_run_bar(self) -> QHBoxLayout:
        bar = QHBoxLayout()
        self._run_label = QLabel("<i>no run attached</i>")
        self._run_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        open_button = QPushButton("Open run…")
        open_button.clicked.connect(self.choose_run)
        self._detach_button = QPushButton("Detach")
        self._detach_button.clicked.connect(self.detach_run)
        self._detach_button.setEnabled(False)

        bar.addWidget(QLabel("<b>Run:</b>"))
        bar.addWidget(self._run_label, stretch=1)
        bar.addWidget(open_button)
        bar.addWidget(self._detach_button)
        return bar

    def _build_progress(self) -> QWidget:
        box = QWidget()
        layout = QVBoxLayout(box)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)

        self._progress = QProgressBar()
        self._progress.setTextVisible(True)
        self._progress.setFormat("%v / %m  (%p%)")

        self._summary = QLabel("—")
        self._summary.setTextFormat(Qt.TextFormat.RichText)

        layout.addWidget(self._progress)
        layout.addWidget(self._summary)
        return box

    def _build_tabs(self) -> QWidget:
        tabs = QTabWidget()

        self._table = QTableView()
        self._table.setModel(self._model)
        self._table.setSortingEnabled(True)
        self._table.setSelectionBehavior(QTableView.SelectionBehavior.SelectRows)
        self._table.setAlternatingRowColors(True)
        self._table.verticalHeader().setVisible(False)
        # Uniform row heights let Qt compute geometry without measuring every
        # row -- required for the view to stay responsive at 47,000 rows.
        self._table.setVerticalScrollMode(QTableView.ScrollMode.ScrollPerPixel)
        self._table.verticalHeader().setDefaultSectionSize(22)
        header = self._table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        header.setStretchLastSection(True)

        # setSortingEnabled() immediately sorts by the header's current indicator,
        # which defaults to column 0 descending -- and since column 0 (rank) maps
        # onto affinity, that silently presents the *worst* binders first. State
        # the intended order explicitly so the header and the model agree: most
        # negative affinity, i.e. the best hits, at the top.
        self._table.sortByColumn(AFFINITY_COLUMN, Qt.SortOrder.AscendingOrder)

        self._ligands = LigandPanel()
        self._setup = SetupPanel()

        # Finishing preparation fills in the ligand directory on the Setup tab and
        # moves the user there: the two stages are separate processes but one
        # workflow, and re-typing the path is the sort of friction that leads to
        # screening the wrong directory.
        self._ligands.prepared.connect(self._on_ligands_prepared)
        # Launching a screen attaches this window to the run it just started, so
        # the user lands on results rather than having to find the directory.
        self._setup.runStarted.connect(self._on_run_started)

        self._tabs = tabs
        # Ordered as the work is done: prepare, then screen, then read results.
        tabs.addTab(_scrolled(self._ligands), "1. Ligands")
        tabs.addTab(_scrolled(self._setup), "2. Screen")
        tabs.addTab(self._table, "3. Results")
        return tabs

    def _on_ligands_prepared(self, ligand_dir: str) -> None:
        self._setup.ligands.set_path(ligand_dir)
        self._tabs.setCurrentIndex(1)
        self._log_line(f"ligands prepared: {ligand_dir}")

    def _on_run_started(self, run_dir: str) -> None:
        self.attach_run(run_dir)
        self._tabs.setCurrentWidget(self._table)

    def _build_activity(self) -> QWidget:
        self._log = QPlainTextEdit()
        self._log.setReadOnly(True)
        # Qt enforces the cap itself, so the widget cannot grow without bound
        # even if something upstream forgets to trim.
        self._log.setMaximumBlockCount(ACTIVITY_LINES)
        self._log.setFont(QFont("monospace", 9))
        self._log.setPlaceholderText("Activity will appear here once a run is attached.")
        return self._log

    def _build_menu(self) -> None:
        file_menu = self.menuBar().addMenu("&File")

        open_action = QAction("&Open run…", self)
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self.choose_run)
        file_menu.addAction(open_action)

        file_menu.addSeparator()
        quit_action = QAction("&Quit", self)
        quit_action.setShortcut("Ctrl+Q")
        quit_action.triggered.connect(self.close)
        file_menu.addAction(quit_action)

    def _connect(self) -> None:
        self._watcher.statusChanged.connect(self._on_status)
        self._watcher.recordsAdded.connect(self._on_records)
        self._watcher.runFinished.connect(self._on_finished)
        self._watcher.error.connect(self._on_error)

    # -- actions -----------------------------------------------------------

    def choose_run(self) -> None:
        directory = QFileDialog.getExistingDirectory(self, "Open run directory")
        if directory:
            self.attach_run(directory)

    def attach_run(self, run_dir: str) -> None:
        self._model.reset([])
        self._activity.clear()
        self._log.clear()

        if not self._watcher.attach(run_dir):
            return

        # Replay an existing journal in one go rather than trickling it through
        # the poll loop, which would take minutes for a large finished run.
        existing = self._watcher.load_all()
        if existing:
            self._model.reset(existing)
            self._log_line(f"loaded {len(existing)} existing records")

        self._run_label.setText(str(Path(run_dir).resolve()))
        self._detach_button.setEnabled(True)
        self.statusBar().showMessage(f"Watching {run_dir}")
        self._log_line(f"attached to {run_dir}")

    def detach_run(self) -> None:
        self._watcher.detach()
        self._run_label.setText("<i>no run attached</i>")
        self._detach_button.setEnabled(False)
        self.statusBar().showMessage("Detached")
        self._log_line("detached")

    # -- watcher signals ---------------------------------------------------

    def _on_status(self, status: RunStatus) -> None:
        self._progress.setMaximum(max(status.total, 1))
        self._progress.setValue(status.done)

        parts = [
            f"<b>{status.state}</b>",
            f"ok <b>{status.completed}</b>",
            f"failed <b>{status.failed}</b>",
        ]
        if status.engine:
            parts.append(f"engine <b>{status.engine}</b>")
        if (rate := status.rate_per_s) is not None:
            parts.append(f"{rate:.2f} lig/s")
        if (eta := status.eta_s) is not None:
            parts.append(f"eta <b>{_format_duration(eta)}</b>")
        self._summary.setText("  ·  ".join(parts))

    def _on_records(self, records: list[LigandResult]) -> None:
        self._model.append(records)
        for record in records[-20:]:
            if record.status == "ok":
                affinity = record.best_affinity
                shown = "n/a" if affinity is None else f"{affinity:.1f}"
                self._log_line(f"{record.ligand_id}  {shown} kcal/mol  ({record.elapsed_s:.1f}s)")
            else:
                self._log_line(f"{record.ligand_id}  {record.status}: {record.error or ''}")
        if len(records) > 20:
            self._log_line(f"… and {len(records) - 20} more")

    def _on_finished(self, status: RunStatus) -> None:
        self._log_line(f"run {status.state}: {status.completed} ok, {status.failed} failed")
        self.statusBar().showMessage(f"Run {status.state}")

    def _on_error(self, message: str) -> None:
        self._log_line(f"error: {message}")
        self.statusBar().showMessage(message)

    # -- helpers -----------------------------------------------------------

    def _log_line(self, text: str) -> None:
        self._activity.append(text)
        self._log.appendPlainText(text)


def _scrolled(widget: QWidget) -> QScrollArea:
    """Wrap a panel so it stays usable on a short window."""
    area = QScrollArea()
    area.setWidget(widget)
    area.setWidgetResizable(True)
    return area


def _format_duration(seconds: float) -> str:
    seconds = int(seconds)
    hours, rem = divmod(seconds, 3600)
    minutes, secs = divmod(rem, 60)
    if hours:
        return f"{hours}h {minutes}m"
    if minutes:
        return f"{minutes}m {secs}s"
    return f"{secs}s"


def main(argv: list[str] | None = None) -> int:
    """Launch the GUI. An optional run directory argument is attached on start."""
    argv = list(sys.argv if argv is None else argv)
    app = QApplication(argv)
    app.setApplicationName("Drydock")

    window = MainWindow()
    window.show()

    if len(argv) > 1 and Path(argv[1]).is_dir():
        # Deferred so the window is painted before any journal reading begins.
        QTimer.singleShot(0, lambda: window.attach_run(argv[1]))

    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
