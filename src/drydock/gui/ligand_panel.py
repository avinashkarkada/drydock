"""The Ligands panel: convert, protonate, optimise and prepare a library.

This is the first half of the tool, and the half that decides what the second
half is actually docking. Protonation state, ring geometry and macrocycle
handling all change results, and none of them announce themselves afterwards --
a PDBQT built from a badly protonated molecule looks exactly like one built well.

So the panel surfaces those choices rather than burying them behind a default,
and explains the one that most often surprises people: Vina samples rotatable
bonds but never ring conformations, so whatever ring geometry is in the input is
held rigid for the entire docking run.

Like the Setup panel, preparation is launched as a detached process and watched
through the directory it writes. Closing the window does not stop a preparation
run that may take half an hour.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from drydock.gui.setup_panel import PathPicker, _drydock_executable

POLL_MS = 1000


class LigandPanel(QWidget):
    """Prepare a compound library for docking."""

    prepared = Signal(str)  # prepared ligand directory

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._process: subprocess.Popen | None = None
        self._expected = 0
        # Guards completion handling. _finish() clears _process, which the
        # completion check would otherwise read as "finished" a second time --
        # re-emitting `prepared` and re-running whatever is connected to it.
        self._running = False

        self._timer = QTimer(self)
        self._timer.setInterval(POLL_MS)
        self._timer.timeout.connect(self._poll)

        self._build()

    def _build(self) -> None:
        layout = QVBoxLayout(self)
        layout.setSpacing(10)

        layout.addWidget(self._input_group())
        layout.addWidget(self._chemistry_group())
        layout.addWidget(self._filters_group())

        self._progress = QProgressBar()
        self._progress.setFormat("%v / %m  (%p%)")
        self._progress.setVisible(False)
        layout.addWidget(self._progress)

        self._status = QLabel()
        self._status.setWordWrap(True)
        self._status.setTextFormat(Qt.TextFormat.RichText)
        layout.addWidget(self._status)

        buttons = QHBoxLayout()
        self._survey_button = QPushButton("Survey library")
        self._survey_button.clicked.connect(self.survey_library)
        self._prepare_button = QPushButton("Prepare ligands")
        self._prepare_button.setDefault(True)
        self._prepare_button.clicked.connect(self.prepare)
        self._stop_button = QPushButton("Stop")
        self._stop_button.clicked.connect(self.stop)
        self._stop_button.setEnabled(False)

        buttons.addStretch(1)
        buttons.addWidget(self._survey_button)
        buttons.addWidget(self._stop_button)
        buttons.addWidget(self._prepare_button)
        layout.addLayout(buttons)
        layout.addStretch(1)

    def _input_group(self) -> QGroupBox:
        group = QGroupBox("Library")
        form = QFormLayout(group)

        self.library = PathPicker("SDF, SMILES or MOL2 (.gz accepted)")
        self.out_dir = PathPicker("Directory to write PDBQTs and the manifest", directory=True)

        self.id_field = QLineEdit()
        self.id_field.setPlaceholderText("auto-detected (COMPOUND_ID, ChEMBL_ID, Name, …)")

        form.addRow("Input file", self.library)
        form.addRow("Output directory", self.out_dir)
        form.addRow("ID field", self.id_field)

        note = QLabel(
            "<b>Survey library</b> first: identifiers repeat in real libraries. "
            "CMNPD 1.0 has 47,451 records under 25,224 identifiers because every "
            "stereoisomer is enumerated under one ID. Drydock keeps them all under "
            "distinct filenames, but the count is worth knowing before you start."
        )
        note.setWordWrap(True)
        note.setTextFormat(Qt.TextFormat.RichText)
        form.addRow(note)
        return group

    def _chemistry_group(self) -> QGroupBox:
        group = QGroupBox("Chemistry")
        form = QFormLayout(group)

        self.ph = QDoubleSpinBox()
        self.ph.setRange(0.0, 14.0)
        self.ph.setSingleStep(0.1)
        self.ph.setValue(7.4)
        self.ph.setDecimals(1)

        self.geometry = QComboBox()
        self.geometry.addItems(
            [
                "Generate 3D and minimise (MMFF94s)",
                "Keep input coordinates",
            ]
        )

        self.conformers = QSpinBox()
        self.conformers.setRange(1, 50)
        self.conformers.setValue(1)

        self.macrocycles = QCheckBox("Let Vina open and re-close macrocyclic rings")
        self.macrocycles.setChecked(True)

        self.tautomers = QCheckBox("Enumerate tautomers (multiplies the library)")

        form.addRow("pH", self.ph)
        form.addRow("Geometry", self.geometry)
        form.addRow("Conformers per compound", self.conformers)
        form.addRow(self.macrocycles)
        form.addRow(self.tautomers)

        note = QLabel(
            "<b>Ring conformations:</b> Vina samples rotatable bonds but never ring "
            "geometry — whatever ring shape is in the PDBQT is held rigid for the "
            "whole run. Minimising repairs bad geometry but only finds the nearest "
            "minimum; it will not flip a chair to a boat. More than one conformer "
            "is the only setting here that genuinely samples ring space, at "
            "proportionally more docking time."
        )
        note.setWordWrap(True)
        note.setTextFormat(Qt.TextFormat.RichText)
        form.addRow(note)
        return group

    def _filters_group(self) -> QGroupBox:
        group = QGroupBox("Filters and execution")
        form = QFormLayout(group)

        self.max_torsions = QSpinBox()
        self.max_torsions.setRange(0, 100)
        self.max_torsions.setValue(0)
        self.max_torsions.setSpecialValueText("no limit")

        self.limit = QSpinBox()
        self.limit.setRange(0, 10_000_000)
        self.limit.setValue(0)
        self.limit.setSpecialValueText("whole library")

        self.workers = QSpinBox()
        self.workers.setRange(0, 256)
        self.workers.setValue(0)
        self.workers.setSpecialValueText("all cores")

        self.seed = QSpinBox()
        self.seed.setRange(0, 2**31 - 1)
        self.seed.setValue(0)

        self.resume = QCheckBox("Skip ligands already prepared")
        self.resume.setChecked(True)

        form.addRow("Max rotatable bonds", self.max_torsions)
        form.addRow("Only prepare first N", self.limit)
        form.addRow("Parallel jobs", self.workers)
        form.addRow("Seed", self.seed)
        form.addRow(self.resume)
        return group

    # -- actions -----------------------------------------------------------

    def survey_library(self) -> None:
        """Report records against distinct compounds, without preparing anything."""
        path = self.library.path()
        if not path:
            self._note("Choose a library file first.", error=True)
            return

        from drydock.core.library import LibraryFormatError, survey

        self._note("Surveying…")
        try:
            info = survey(path, id_field=self.id_field.text().strip() or None)
        except (LibraryFormatError, OSError) as exc:
            self._note(f"Could not read the library: {exc}", error=True)
            return

        self._expected = info["records"]
        lines = [
            f"<b>{info['records']:,}</b> records, "
            f"<b>{info['distinct_compounds']:,}</b> distinct compounds",
        ]
        if info["compounds_with_variants"]:
            lines.append(
                f"{info['compounds_with_variants']:,} identifiers repeat "
                f"(up to {info['max_variants']} times) — these are kept as separate "
                "ligands and grouped again in the results."
            )
        self._note("<br>".join(lines))

    def prepare(self) -> None:
        """Launch preparation as a detached process."""
        command = self._build_command()
        if command is None:
            return

        out = Path(self.out_dir.path())
        try:
            out.mkdir(parents=True, exist_ok=True)
            log = out / "prep.log"
            with open(log, "ab") as handle:
                self._process = subprocess.Popen(  # noqa: S603 - argv is constructed
                    command,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    stdin=subprocess.DEVNULL,
                    start_new_session=True,
                    cwd=os.getcwd(),
                )
        except OSError as exc:
            QMessageBox.critical(self, "Could not start", f"Failed to launch:\n{exc}")
            return

        if not self._expected:
            # Cheap, and turns an indeterminate bar into a real one.
            try:
                from drydock.core.library import count_records

                self._expected = count_records(self.library.path())
            except Exception:  # noqa: BLE001 - a progress total is not worth failing over
                self._expected = 0

        self._running = True
        self._progress.setVisible(True)
        self._progress.setMaximum(max(self._expected, 1))
        self._progress.setValue(0)
        self._prepare_button.setEnabled(False)
        self._stop_button.setEnabled(True)
        self._note(f"Preparing into <tt>{out}</tt>. This window can be closed.")
        self._timer.start()

    def stop(self) -> None:
        """Ask the preparation process to stop.

        Preparation is resumable, so stopping costs only the ligands in flight;
        restarting with the same output directory picks up where it left off.
        """
        if self._process and self._process.poll() is None:
            self._process.terminate()
            self._note("Stopping… already-prepared ligands are kept and can be resumed.")
        self._finish()

    def _build_command(self) -> list[str] | None:
        library, out = self.library.path(), self.out_dir.path()
        if not library:
            self._note("Choose a library file.", error=True)
            return None
        if not out:
            self._note("Choose an output directory.", error=True)
            return None

        command = [
            _drydock_executable(),
            "prep-ligands",
            "--input", library,
            "--out", out,
            "--ph", f"{self.ph.value():.1f}",
            "--conformers", str(self.conformers.value()),
            "--workers", str(self.workers.value()),
            "--seed", str(self.seed.value()),
        ]
        command.append(
            "--optimize" if self.geometry.currentIndex() == 0 else "--no-optimize"
        )
        command.append(
            "--macrocycles" if self.macrocycles.isChecked() else "--rigid-macrocycles"
        )
        if self.tautomers.isChecked():
            command.append("--tautomers")
        if self.max_torsions.value():
            command += ["--max-torsions", str(self.max_torsions.value())]
        if self.limit.value():
            command += ["--limit", str(self.limit.value())]
        if self.id_field.text().strip():
            command += ["--id-field", self.id_field.text().strip()]
        if not self.resume.isChecked():
            command.append("--no-resume")
        return command

    # -- polling -----------------------------------------------------------

    def _poll(self) -> None:
        """Read progress from the output directory, not from the process."""
        if not self._running:
            return

        from drydock.core.prep_runner import PrepDir

        prep = PrepDir(self.out_dir.path())
        prepared, failed = prep.progress()

        self._progress.setMaximum(max(self._expected, prepared + failed, 1))
        self._progress.setValue(prepared + failed)

        message = f"prepared <b>{prepared:,}</b>"
        if failed:
            message += f", failed <b>{failed:,}</b>"
        self._note(message)

        finished = self._process is None or self._process.poll() is not None
        if finished:
            info = prep.read_info()
            summary = (info or {}).get("summary", {})
            detail = (
                f"prepared <b>{summary.get('prepared', prepared):,}</b>, "
                f"failed <b>{summary.get('failed', failed):,}</b>"
            )
            if summary.get("failed"):
                detail += f" — see <tt>{prep.failures_file}</tt>"
            self._note(f"Finished: {detail}<br>Manifest: <tt>{prep.manifest_file}</tt>")
            self._finish()
            self.prepared.emit(str(prep.path))

    def _finish(self) -> None:
        self._running = False
        self._timer.stop()
        self._process = None
        self._prepare_button.setEnabled(True)
        self._stop_button.setEnabled(False)

    def _note(self, html: str, error: bool = False) -> None:
        colour = "#c0392b" if error else "#555"
        self._status.setText(f"<div style='color:{colour}'>{html}</div>")
