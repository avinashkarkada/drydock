"""The Receptor panel: prepare a structure, and check what came out.

Preparation is fast enough (seconds) to run in-process rather than detached, so
this panel is simpler than the other two. What it does not do is report success
and stop there: the result is inspected and the findings shown, because a
receptor can prepare cleanly and still be unusable.

The case that motivated it: a receptor prepared elsewhere reached Drydock with
zero polar hydrogens. Every hydrogen-bond donor in the protein was missing --
backbone amides, lysines, tyrosine hydroxyls, all scoring acceptor-only. Nothing
about the file looked wrong, and no error was raised at any point.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from drydock.gui.setup_panel import PathPicker


class ReceptorPanel(QWidget):
    """Prepare a receptor structure and report on the result."""

    prepared = Signal(str)  # prepared receptor PDBQT

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._build()

    def _build(self) -> None:
        layout = QVBoxLayout(self)
        layout.setSpacing(10)

        layout.addWidget(self._input_group())
        layout.addWidget(self._options_group())
        layout.addWidget(self._zinc_group())

        self._status = QLabel()
        self._status.setWordWrap(True)
        self._status.setTextFormat(Qt.TextFormat.RichText)
        self._status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(self._status)

        buttons = QHBoxLayout()
        self._check_button = QPushButton("Check existing PDBQT")
        self._check_button.clicked.connect(self.check_existing)
        self._prepare_button = QPushButton("Prepare receptor")
        self._prepare_button.setDefault(True)
        self._prepare_button.clicked.connect(self.prepare)

        buttons.addStretch(1)
        buttons.addWidget(self._check_button)
        buttons.addWidget(self._prepare_button)
        layout.addLayout(buttons)
        layout.addStretch(1)

    def _input_group(self) -> QGroupBox:
        group = QGroupBox("Structure")
        form = QFormLayout(group)

        self.structure = PathPicker("PDB or mmCIF (.pdb, .cif)")
        # A save target, not an existing file, and a basename rather than a
        # filename -- so the dialog must allow naming something new, and any
        # .pdbqt the user supplies is dropped rather than doubled.
        self.out = PathPicker(
            "Output name; .pdbqt is added for you",
            save=True,
            strip_suffix=".pdbqt",
        )

        form.addRow("Input structure", self.structure)
        form.addRow("Output", self.out)

        note = QLabel(
            "Adds hydrogens, assigns AutoDock atom types and partial charges. "
            "mmCIF is converted automatically. The result is checked before it is "
            "handed on — a receptor can prepare cleanly and still be unusable."
        )
        note.setWordWrap(True)
        form.addRow(note)
        return group

    def _options_group(self) -> QGroupBox:
        group = QGroupBox("Options")
        form = QFormLayout(group)

        self.allow_bad = QCheckBox("Delete residues with missing atoms")
        self.allow_bad.setChecked(True)

        self.delete_residues = QLineEdit()
        self.delete_residues.setPlaceholderText("e.g. A:350,B:15,16 — waters, artefacts")

        self.charge_model = QComboBox()
        self.charge_model.addItems(["gasteiger", "espaloma", "zero"])

        self.altloc = QLineEdit()
        self.altloc.setPlaceholderText("e.g. A — only if the structure has alternates")
        self.altloc.setMaximumWidth(120)

        form.addRow(self.allow_bad)
        form.addRow("Delete residues", self.delete_residues)
        form.addRow("Charges", self.charge_model)
        form.addRow("Alternate location", self.altloc)

        note = QLabel(
            "Crystal structures routinely have disordered side chains, so deleting "
            "incomplete residues is on by default — refusing a whole structure over "
            "one unresolved lysine helps nobody. Anything removed is listed below."
        )
        note.setWordWrap(True)
        form.addRow(note)
        return group

    def _zinc_group(self) -> QGroupBox:
        group = QGroupBox("Zinc metalloproteins (AutoDock4Zn)")
        layout = QVBoxLayout(group)

        self.add_zinc = QCheckBox(
            "Add tetrahedral zinc pseudo-atoms after preparing (needed for the ad4 engine)"
        )
        layout.addWidget(self.add_zinc)

        note = QLabel(
            "Pseudo-atoms go only at <i>vacant</i> coordination sites. A zinc already "
            "saturated by the protein — a structural rather than catalytic site — "
            "correctly gets none, since there is nowhere for a ligand to bind."
        )
        note.setWordWrap(True)
        note.setTextFormat(Qt.TextFormat.RichText)
        layout.addWidget(note)
        return group

    # -- actions -----------------------------------------------------------

    def check_existing(self) -> None:
        """Inspect a PDBQT without preparing anything."""
        path = self.structure.path() or self.out.path()
        if not path:
            self._note("Choose a file to check.", error=True)
            return

        from drydock.core.receptor import inspect

        self._show_report(inspect(path), Path(path))

    def prepare(self) -> None:
        structure = self.structure.path()
        out = self.out.path()
        if not structure:
            self._note("Choose an input structure.", error=True)
            return
        if not out:
            self._note("Choose an output path.", error=True)
            return

        from drydock.core.recprep import ReceptorPrepError, prepare_receptor

        self._prepare_button.setEnabled(False)
        self._note("Preparing…")
        try:
            result = prepare_receptor(
                structure,
                out,
                allow_bad_residues=self.allow_bad.isChecked(),
                delete_residues=self.delete_residues.text().strip() or None,
                charge_model=self.charge_model.currentText(),
                default_altloc=self.altloc.text().strip() or None,
            )
        except ReceptorPrepError as exc:
            self._note(str(exc).replace("\n", "<br>"), error=True)
            return
        finally:
            self._prepare_button.setEnabled(True)

        receptor = result.receptor_pdbqt
        extra = list(result.deleted_residues)

        if self.add_zinc.isChecked():
            from drydock.core.zinc import ZincError, add_zinc_pseudo_atoms

            try:
                zinc = add_zinc_pseudo_atoms(receptor)
            except ZincError as exc:
                extra.append(f"zinc pseudo-atoms not added: {exc}")
            else:
                receptor = zinc.receptor_tz
                extra.append(
                    f"placed {zinc.n_pseudo_atoms} zinc pseudo-atom(s) "
                    f"around {zinc.n_zinc} zinc(s)"
                )

        from drydock.core.receptor import inspect

        self._show_report(inspect(receptor), receptor, extra=extra)
        self.prepared.emit(str(receptor))

    def _show_report(self, report, path: Path, extra: list[str] | None = None) -> None:
        lines = [
            f"<b>{path}</b>",
            f"{report.n_atoms} atoms &middot; "
            f"chains {', '.join(report.chains) or '-'} &middot; "
            f"polar H <b>{report.n_polar_hydrogens}</b> &middot; "
            f"metals {', '.join(report.metals) or 'none'}",
        ]
        lines += [f"<i>{item}</i>" for item in (extra or [])]
        lines += [f"note: {note}" for note in report.notes]
        lines += [
            f"<b style='color:#c0392b'>problem: {problem}</b>" for problem in report.problems
        ]
        if report.ok:
            lines.append("<b style='color:#207d3f'>receptor looks usable</b>")

        self._note("<br>".join(lines), error=bool(report.problems))

    def _note(self, html: str, error: bool = False) -> None:
        colour = "#c0392b" if error else "#555"
        self._status.setText(f"<div style='color:{colour}'>{html}</div>")
