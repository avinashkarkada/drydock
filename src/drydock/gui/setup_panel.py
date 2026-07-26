"""The Setup panel: configure a screen and launch it.

Launching is the one place the GUI could compromise the design, so it does not
run anything in-process. It writes a configuration, spawns ``drydock screen`` as a
**detached** process, and then attaches to the resulting run directory as an
ordinary watcher.

The consequence is the property the whole architecture exists for: closing this
window does not stop the run, and a run that dies does not take the window with
it. The GUI has no privileged relationship with the work -- it is simply the first
reader of a run directory that anything else could equally read.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)


class PathPicker(QWidget):
    """A line edit with a browse button."""

    changed = Signal(str)

    def __init__(self, placeholder: str = "", directory: bool = False, parent=None) -> None:
        super().__init__(parent)
        self._directory = directory

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.edit = QLineEdit()
        self.edit.setPlaceholderText(placeholder)
        self.edit.textChanged.connect(self.changed.emit)

        button = QPushButton("Browse…")
        button.clicked.connect(self._browse)

        layout.addWidget(self.edit, stretch=1)
        layout.addWidget(button)

    def _browse(self) -> None:
        if self._directory:
            chosen = QFileDialog.getExistingDirectory(self, "Select directory")
        else:
            chosen, _ = QFileDialog.getOpenFileName(self, "Select file")
        if chosen:
            self.edit.setText(chosen)

    def path(self) -> str:
        return self.edit.text().strip()

    def set_path(self, value: str) -> None:
        self.edit.setText(value)


class SetupPanel(QWidget):
    """Configure and launch a screen."""

    runStarted = Signal(str)  # run directory

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._build()

    def _build(self) -> None:
        layout = QVBoxLayout(self)
        layout.setSpacing(10)

        layout.addWidget(self._inputs_group())
        layout.addWidget(self._box_group())
        layout.addWidget(self._engine_group())

        self._status = QLabel()
        self._status.setWordWrap(True)
        self._status.setTextFormat(Qt.TextFormat.RichText)
        layout.addWidget(self._status)

        buttons = QHBoxLayout()
        self._check_button = QPushButton("Check receptor")
        self._check_button.clicked.connect(self.check_receptor)
        self._start_button = QPushButton("Start screen")
        self._start_button.setDefault(True)
        self._start_button.clicked.connect(self.start_screen)

        buttons.addStretch(1)
        buttons.addWidget(self._check_button)
        buttons.addWidget(self._start_button)
        layout.addLayout(buttons)
        layout.addStretch(1)

    def _inputs_group(self) -> QGroupBox:
        group = QGroupBox("Inputs")
        form = QFormLayout(group)

        self.receptor = PathPicker("Prepared receptor PDBQT")
        self.ligands = PathPicker("Prepared ligand directory", directory=True)
        self.run_dir = PathPicker("Run directory to create", directory=True)

        form.addRow("Receptor", self.receptor)
        form.addRow("Ligands", self.ligands)
        form.addRow("Run directory", self.run_dir)

        note = QLabel(
            "Receptors are not prepared here -- point Drydock at a PDBQT you "
            "prepared yourself. Use <b>Check receptor</b> first: a receptor missing "
            "polar hydrogens has no hydrogen-bond donors at all, and nothing else "
            "will report that."
        )
        note.setWordWrap(True)
        note.setTextFormat(Qt.TextFormat.RichText)
        form.addRow(note)
        return group

    def _box_group(self) -> QGroupBox:
        group = QGroupBox("Search box")
        form = QFormLayout(group)

        self.residues = QLineEdit()
        self.residues.setPlaceholderText("e.g. 187,188,397,401,405,411,420,421,423")

        self.chain = QLineEdit()
        self.chain.setPlaceholderText("optional, e.g. B")
        self.chain.setMaximumWidth(80)

        self.padding = QSpinBox()
        self.padding.setRange(0, 30)
        self.padding.setValue(5)
        self.padding.setSuffix(" Å")

        self.center = QLineEdit()
        self.center.setPlaceholderText("x,y,z — leave blank to use residues")
        self.size = QLineEdit()
        self.size.setPlaceholderText("x,y,z — leave blank to use residues")

        form.addRow("Active-site residues", self.residues)
        form.addRow("Chain", self.chain)
        form.addRow("Padding", self.padding)
        form.addRow("Or explicit centre", self.center)
        form.addRow("Or explicit size", self.size)
        return group

    def _engine_group(self) -> QGroupBox:
        group = QGroupBox("Docking")
        form = QFormLayout(group)

        self.engine = QComboBox()
        self.engine.addItems(["vina", "vinardo", "ad4"])
        self.engine.currentTextChanged.connect(self._engine_changed)

        self.maps = PathPicker("AutoGrid maps directory (ad4 only)", directory=True)
        self.maps.setEnabled(False)

        self.exhaustiveness = QSpinBox()
        self.exhaustiveness.setRange(1, 512)
        self.exhaustiveness.setValue(8)

        self.modes = QSpinBox()
        self.modes.setRange(1, 50)
        self.modes.setValue(9)

        self.seed = QSpinBox()
        self.seed.setRange(0, 2**31 - 1)
        self.seed.setValue(0)

        self.workers = QSpinBox()
        self.workers.setRange(0, 256)
        self.workers.setValue(0)
        self.workers.setSpecialValueText("all cores")

        self.resume = QCheckBox("Resume if the run directory already has results")
        self.resume.setChecked(True)

        form.addRow("Engine", self.engine)
        form.addRow("Maps (ad4)", self.maps)
        form.addRow("Exhaustiveness", self.exhaustiveness)
        form.addRow("Binding modes", self.modes)
        form.addRow("Seed", self.seed)
        form.addRow("Parallel jobs", self.workers)
        form.addRow(self.resume)
        return group

    def _engine_changed(self, engine: str) -> None:
        needs_maps = engine == "ad4"
        self.maps.setEnabled(needs_maps)
        if needs_maps:
            self._note(
                "The <b>ad4</b> engine uses AutoDock4/AutoDock4Zn scoring over "
                "pre-computed AutoGrid maps. For a zinc metalloprotein, run "
                "<tt>drydock add-zinc-pseudo</tt> then <tt>drydock maps</tt> first."
            )

    # -- actions -----------------------------------------------------------

    def check_receptor(self) -> None:
        receptor = self.receptor.path()
        if not receptor:
            self._note("Choose a receptor first.", error=True)
            return

        from drydock.core.receptor import inspect

        report = inspect(receptor)
        lines = [
            f"{report.n_atoms} atoms, chains {', '.join(report.chains) or '-'}",
            f"polar hydrogens: {report.n_polar_hydrogens}",
            f"metals: {', '.join(report.metals) or 'none'}",
        ]
        for note in report.notes:
            lines.append(f"note: {note}")
        for problem in report.problems:
            lines.append(f"<b style='color:#c0392b'>problem: {problem}</b>")

        self._note("<br>".join(lines), error=bool(report.problems))

    def start_screen(self) -> None:
        """Launch a detached screening process and hand its directory back."""
        command = self._build_command()
        if command is None:
            return

        try:
            self._spawn_detached(command)
        except OSError as exc:
            QMessageBox.critical(self, "Could not start", f"Failed to launch:\n{exc}")
            return

        run_dir = self.run_dir.path()
        self._note(
            f"Started. Watching <tt>{run_dir}</tt>.<br>"
            "This window can be closed at any time -- the run continues, and "
            "reopening it will pick the run back up."
        )
        self.runStarted.emit(run_dir)

    def _build_command(self) -> list[str] | None:
        receptor, ligands, run_dir = (
            self.receptor.path(),
            self.ligands.path(),
            self.run_dir.path(),
        )

        missing = [
            name
            for name, value in (
                ("receptor", receptor),
                ("ligand directory", ligands),
                ("run directory", run_dir),
            )
            if not value
        ]
        if missing:
            self._note(f"Still needed: {', '.join(missing)}.", error=True)
            return None

        has_residues = bool(self.residues.text().strip())
        has_explicit = bool(self.center.text().strip() and self.size.text().strip())
        if not (has_residues or has_explicit):
            self._note(
                "Define the box: either active-site residues, or both an explicit "
                "centre and size.",
                error=True,
            )
            return None

        if self.engine.currentText() == "ad4" and not self.maps.path():
            self._note("The ad4 engine needs a maps directory.", error=True)
            return None

        command = [
            _drydock_executable(),
            "screen",
            "--receptor", receptor,
            "--ligands", ligands,
            "--run", run_dir,
            "--engine", self.engine.currentText(),
            "--exhaustiveness", str(self.exhaustiveness.value()),
            "--modes", str(self.modes.value()),
            "--seed", str(self.seed.value()),
            "--workers", str(self.workers.value()),
        ]
        if has_explicit:
            command += ["--center", self.center.text().strip(),
                        "--size", self.size.text().strip()]
        else:
            command += ["--residues", self.residues.text().strip(),
                        "--pad", str(self.padding.value())]
            if self.chain.text().strip():
                command += ["--chain", self.chain.text().strip()]
        if self.maps.path():
            command += ["--maps", self.maps.path()]
        if not self.resume.isChecked():
            command.append("--no-resume")
        return command

    @staticmethod
    def _spawn_detached(command: list[str]) -> None:
        """Start a process that outlives this one.

        ``start_new_session`` puts the child in its own session, so it does not
        share the GUI's process group and is not signalled when the GUI exits or
        the terminal that launched it closes. Output goes to the run's own log
        rather than to inherited descriptors -- inheriting stdout is how twelve
        docking workers previously ended up writing over each other.
        """
        run_dir = Path(command[command.index("--run") + 1])
        run_dir.mkdir(parents=True, exist_ok=True)
        log = run_dir / "screen.log"

        with open(log, "ab") as handle:
            subprocess.Popen(  # noqa: S603 - argument list is constructed, not shell
                command,
                stdout=handle,
                stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL,
                start_new_session=True,
                cwd=os.getcwd(),
            )

    def _note(self, html: str, error: bool = False) -> None:
        colour = "#c0392b" if error else "#555"
        self._status.setText(f"<div style='color:{colour}'>{html}</div>")


def _drydock_executable() -> str:
    """Locate the drydock CLI belonging to this interpreter.

    Resolved from sys.executable rather than from PATH so the GUI cannot launch a
    different installation than the one it is running from -- which would quietly
    produce results from a different set of pinned dependencies.
    """
    candidate = Path(sys.executable).parent / "drydock"
    if candidate.exists():
        return str(candidate)
    return "drydock"
