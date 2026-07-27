"""The Setup panel: configure a screen and launch it.

Launching is the one place the GUI could compromise the design, so it does not
run anything in-process. It writes a configuration, spawns ``drydock screen`` as a
**detached** process, and then attaches to the resulting run directory as an
ordinary watcher.

The consequence is the property the whole architecture exists for: closing this
window does not stop the run, and a run that dies does not take the window with
it. The GUI has no privileged relationship with the work, it is simply the first
reader of a run directory that anything else could equally read.
"""

from __future__ import annotations

import os
import shutil
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
    """A line edit with a browse button.

    Three modes, because a file dialog that offers the wrong operation is worse
    than no dialog: an *open* dialog cannot name a file that does not exist yet,
    so using one for an output leaves the user unable to do the obvious thing.

    ``directory``  choose an existing folder
    ``save``       name an output, which may not exist yet
    (neither)      choose an existing file
    """

    changed = Signal(str)

    def __init__(
        self,
        placeholder: str = "",
        directory: bool = False,
        save: bool = False,
        strip_suffix: str | None = None,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self._directory = directory
        self._save = save
        self._strip_suffix = strip_suffix

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.edit = QLineEdit()
        self.edit.setPlaceholderText(placeholder)
        self.edit.textChanged.connect(self.changed.emit)

        button = QPushButton("Browse...")
        button.clicked.connect(self._browse)

        layout.addWidget(self.edit, stretch=1)
        layout.addWidget(button)

    def _browse(self) -> None:
        if self._directory:
            chosen = QFileDialog.getExistingDirectory(self, "Select directory")
        elif self._save:
            chosen, _ = QFileDialog.getSaveFileName(self, "Save as")
        else:
            chosen, _ = QFileDialog.getOpenFileName(self, "Select file")
        if chosen:
            self.set_path(chosen)

    def path(self) -> str:
        return self.edit.text().strip()

    def set_path(self, value: str) -> None:
        """Set the path, removing an extension the caller will add itself.

        Fields that hold a *basename* would otherwise accumulate it twice: a save
        dialog naturally produces "receptor.pdbqt", the tool appends ".pdbqt", and
        preparation writes "receptor.pdbqt.pdbqt" and then fails looking for the
        name it expected. Easy to do by hand as well as through the dialog.
        """
        value = value.strip()
        if self._strip_suffix and value.lower().endswith(self._strip_suffix.lower()):
            value = value[: -len(self._strip_suffix)]
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
            "Prepare these on the <b>Receptor</b> and <b>Ligands</b> tabs, or point "
            "at files you prepared elsewhere. Either way use <b>Check receptor</b> "
            "first: a receptor missing polar hydrogens has no hydrogen-bond donors "
            "at all, and nothing else will report that."
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
        self.center.setPlaceholderText("x,y,z, leave blank to use residues")
        self.size = QLineEdit()
        self.size.setPlaceholderText("x,y,z, leave blank to use residues")

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
        self.maps.changed.connect(self._maps_changed)

        self._maps_button = QPushButton("Generate maps...")
        self._maps_button.setEnabled(False)
        self._maps_button.clicked.connect(self.generate_maps)

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

        maps_row = QHBoxLayout()
        maps_row.setContentsMargins(0, 0, 0, 0)
        maps_row.addWidget(self.maps, stretch=1)
        maps_row.addWidget(self._maps_button)
        maps_widget = QWidget()
        maps_widget.setLayout(maps_row)

        form.addRow("Engine", self.engine)
        form.addRow("Maps (ad4)", maps_widget)
        form.addRow("Exhaustiveness", self.exhaustiveness)
        form.addRow("Binding modes", self.modes)
        form.addRow("Seed", self.seed)
        form.addRow("Parallel jobs", self.workers)
        form.addRow(self.resume)
        return group

    def _engine_changed(self, engine: str) -> None:
        needs_maps = engine == "ad4"
        self.maps.setEnabled(needs_maps)
        self._maps_button.setEnabled(needs_maps)
        if needs_maps:
            self._note(
                "The <b>ad4</b> engine scores against pre-computed AutoGrid maps, "
                "which are separate from the receptor. Choose a directory and press "
                "<b>Generate maps</b>, or point at one you made earlier with "
                "<tt>drydock maps</tt>."
            )
        else:
            self._note("")

    def _maps_changed(self, path: str) -> None:
        """Report on a chosen maps directory as soon as it is picked.

        Checking here rather than at launch matters: pointing at the directory
        holding the receptor is a natural mistake, and left to Vina it surfaces as
        a failure on every ligand in the run rather than once, now.
        """
        if self.engine.currentText() != "ad4" or not path:
            return

        from drydock.core.zinc import maps_status

        ok, detail = maps_status(path)
        if ok:
            self._note(f"Maps look usable, {detail}.")
        else:
            self._note(detail, error=True)

    def generate_maps(self) -> None:
        """Compute AutoGrid maps for the current receptor and box."""
        receptor = self.receptor.path()
        if not receptor:
            self._note("Choose a receptor first.", error=True)
            return

        box = self._resolve_box()
        if box is None:
            return

        target = self.maps.path()
        if not target:
            self._note(
                "Choose a directory for the maps first, somewhere separate from "
                "the receptor, since roughly 50 MB of grid files will be written.",
                error=True,
            )
            return

        from drydock.core.receptor import inspect
        from drydock.core.zinc import ZincError, run_autogrid, write_gpf

        report = inspect(receptor)
        if report.metals and not report.has_zinc_pseudo_atoms:
            self._note(
                "This receptor has metals but no zinc pseudo-atoms, so the maps "
                "will not carry AutoDock4Zn's coordination geometry. Prepare it "
                "again on the Receptor tab with pseudo-atoms enabled.",
                error=True,
            )
            return

        out = Path(target)
        out.mkdir(parents=True, exist_ok=True)
        local = out / Path(receptor).name
        if Path(receptor).resolve() != local.resolve():
            shutil.copy(receptor, local)

        self._note(f"Running autogrid4 for {box}... this takes a few seconds.")
        self._maps_button.setEnabled(False)
        try:
            gpf = write_gpf(local, box, out / "receptor.gpf")
            maps_dir = run_autogrid(gpf, out)
        except ZincError as exc:
            self._note(str(exc).replace("\n", "<br>"), error=True)
            return
        finally:
            self._maps_button.setEnabled(True)

        written = sorted(maps_dir.glob("*.map"))
        size_mb = sum(m.stat().st_size for m in written) / 1e6
        self._note(
            f"Wrote <b>{len(written)}</b> maps ({size_mb:.0f} MB) to <tt>{maps_dir}</tt>. "
            "Ready to screen with the ad4 engine."
        )

    def _resolve_box(self):
        """Build the box from whichever definition the user filled in."""
        from drydock.core.box import Box
        from drydock.core.receptor import box_from_residues

        centre, size = self.center.text().strip(), self.size.text().strip()
        if centre and size:
            try:
                return Box(
                    tuple(float(v) for v in centre.split(",")),
                    tuple(float(v) for v in size.split(",")),
                )
            except ValueError:
                self._note("Centre and size must each be three numbers, x,y,z.", error=True)
                return None

        residues = self.residues.text().strip()
        if not residues:
            self._note(
                "Define the box first: either active-site residues, or both an "
                "explicit centre and size.",
                error=True,
            )
            return None

        try:
            numbers = [int(v) for v in residues.replace(" ", "").split(",") if v]
            box, _ = box_from_residues(
                self.receptor.path(),
                numbers,
                padding=float(self.padding.value()),
                chain=self.chain.text().strip() or None,
            )
        except ValueError as exc:
            self._note(str(exc), error=True)
            return None
        return box

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
            "This window can be closed at any time, the run continues, and "
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

        if self.engine.currentText() == "ad4":
            # Refuse to start rather than let every ligand fail identically. A
            # screen launched against a directory with no maps produces one
            # failure per compound and no usable diagnosis.
            from drydock.core.zinc import maps_status

            ok, detail = maps_status(self.maps.path())
            if not ok:
                self._note(detail, error=True)
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
        rather than to inherited descriptors, inheriting stdout is how twelve
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
    different installation than the one it is running from, which would quietly
    produce results from a different set of pinned dependencies.
    """
    candidate = Path(sys.executable).parent / "drydock"
    if candidate.exists():
        return str(candidate)
    return "drydock"
