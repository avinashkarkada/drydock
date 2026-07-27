"""Tests for the Ligands and Screen panels.

Focused on the command each panel builds and on the completion handling, rather
than on layout. The command is where a wrong flag silently changes the chemistry
of a whole library; the completion handling is where a panel can fire its
"finished" signal more than once, which it did.
"""

from __future__ import annotations

import pytest

from drydock.gui.ligand_panel import LigandPanel
from drydock.gui.setup_panel import SetupPanel


def _flag_value(command: list[str], flag: str) -> str | None:
    return command[command.index(flag) + 1] if flag in command else None


def _fake_maps_dir(path):
    """A directory shaped like one autogrid4 produces."""
    path.mkdir(parents=True, exist_ok=True)
    (path / "receptor.maps.fld").write_text("# fld\n", encoding="utf-8")
    for atom_type in ("C", "A", "N", "OA", "HD"):
        (path / f"receptor.{atom_type}.map").write_text("# map\n", encoding="utf-8")
    return path


class TestMapsValidation:
    """Checking a maps directory before a run rather than during it."""

    def test_accepts_a_complete_map_set(self, tmp_path):
        from drydock.core.zinc import maps_status

        ok, detail = maps_status(_fake_maps_dir(tmp_path / "maps"))
        assert ok
        assert "5 maps" in detail

    def test_rejects_a_receptor_directory_and_says_why(self, tmp_path):
        from drydock.core.zinc import maps_status

        receptors = tmp_path / "Receptor"
        receptors.mkdir()
        (receptors / "2OVX_TZ.pdbqt").write_text("ATOM\n", encoding="utf-8")

        ok, detail = maps_status(receptors)
        assert not ok
        assert "receptor files" in detail
        assert "drydock maps" in detail

    def test_rejects_an_empty_directory(self, tmp_path):
        from drydock.core.zinc import maps_status

        empty = tmp_path / "empty"
        empty.mkdir()
        ok, detail = maps_status(empty)

        assert not ok
        assert "no AutoGrid maps" in detail

    def test_rejects_maps_missing_the_field_file(self, tmp_path):
        """An incomplete set is worse than none: it looks plausible."""
        from drydock.core.zinc import maps_status

        partial = tmp_path / "partial"
        partial.mkdir()
        (partial / "receptor.C.map").write_text("# map\n", encoding="utf-8")

        ok, detail = maps_status(partial)
        assert not ok
        assert "incomplete" in detail

    def test_rejects_nothing_at_all(self):
        from drydock.core.zinc import maps_status

        ok, detail = maps_status(None)
        assert not ok
        assert "ad4" in detail

    def test_rejects_a_path_that_is_not_a_directory(self, tmp_path):
        from drydock.core.zinc import maps_status

        file = tmp_path / "a.pdbqt"
        file.write_text("", encoding="utf-8")

        ok, detail = maps_status(file)
        assert not ok
        assert "not a directory" in detail


class TestLigandPanelCommand:
    @pytest.fixture
    def panel(self, qapp, tmp_path):
        panel = LigandPanel()
        library = tmp_path / "lib.smi"
        library.write_text("CCO ethanol\n", encoding="utf-8")
        panel.library.set_path(str(library))
        panel.out_dir.set_path(str(tmp_path / "out"))
        return panel

    def test_builds_a_command_with_the_defaults(self, panel):
        command = panel._build_command()

        assert command is not None
        assert "prep-ligands" in command
        assert _flag_value(command, "--ph") == "7.4"
        assert "--optimize" in command
        assert "--macrocycles" in command

    def test_missing_library_refuses_and_explains(self, qapp):
        panel = LigandPanel()
        panel.out_dir.set_path("/tmp/out")
        assert panel._build_command() is None
        assert "library" in panel._status.text().lower()

    def test_missing_output_refuses_and_explains(self, qapp, tmp_path):
        panel = LigandPanel()
        panel.library.set_path(str(tmp_path / "lib.smi"))
        assert panel._build_command() is None
        assert "output" in panel._status.text().lower()

    def test_geometry_choice_maps_to_the_right_flag(self, panel):
        """Selecting 'keep input coordinates' must not still regenerate them."""
        panel.geometry.setCurrentIndex(1)
        command = panel._build_command()

        assert "--no-optimize" in command
        assert "--optimize" not in command

    def test_rigid_macrocycles_flag(self, panel):
        panel.macrocycles.setChecked(False)
        command = panel._build_command()

        assert "--rigid-macrocycles" in command
        assert "--macrocycles" not in command

    def test_ph_is_passed_through(self, panel):
        panel.ph.setValue(5.0)
        assert _flag_value(panel._build_command(), "--ph") == "5.0"

    def test_conformers_are_passed_through(self, panel):
        panel.conformers.setValue(4)
        assert _flag_value(panel._build_command(), "--conformers") == "4"

    def test_optional_flags_are_omitted_at_their_defaults(self, panel):
        """A zero-valued spinbox means 'unset', not 'zero'."""
        command = panel._build_command()

        assert "--max-torsions" not in command
        assert "--limit" not in command
        assert "--id-field" not in command
        assert "--tautomers" not in command

    def test_optional_flags_appear_when_set(self, panel):
        panel.max_torsions.setValue(12)
        panel.limit.setValue(500)
        panel.id_field.setText("CATALOG_ID")
        panel.tautomers.setChecked(True)
        command = panel._build_command()

        assert _flag_value(command, "--max-torsions") == "12"
        assert _flag_value(command, "--limit") == "500"
        assert _flag_value(command, "--id-field") == "CATALOG_ID"
        assert "--tautomers" in command

    def test_resume_is_on_by_default_and_can_be_disabled(self, panel):
        assert "--no-resume" not in panel._build_command()

        panel.resume.setChecked(False)
        assert "--no-resume" in panel._build_command()


class TestLigandPanelCompletion:
    def test_prepared_is_not_emitted_before_a_run(self, qapp, tmp_path):
        panel = LigandPanel()
        panel.out_dir.set_path(str(tmp_path))

        emitted: list[str] = []
        panel.prepared.connect(emitted.append)
        panel._poll()

        assert emitted == [], "polling while idle must do nothing"

    def test_prepared_is_emitted_exactly_once(self, qapp, tmp_path):
        """The regression: _finish() clears _process, which the completion check
        then read as 'finished' again, re-emitting on every subsequent poll."""
        out = tmp_path / "out"
        out.mkdir()
        (out / "manifest.csv").write_text("ligand_id\nA\n", encoding="utf-8")
        (out / "prep.json").write_text('{"summary": {"prepared": 1, "failed": 0}}', "utf-8")

        panel = LigandPanel()
        panel.out_dir.set_path(str(out))
        panel._running = True
        panel._process = None

        emitted: list[str] = []
        panel.prepared.connect(emitted.append)

        panel._poll()
        panel._poll()
        panel._poll()

        assert len(emitted) == 1

    def test_finishing_re_enables_the_prepare_button(self, qapp, tmp_path):
        out = tmp_path / "out"
        out.mkdir()
        (out / "manifest.csv").write_text("ligand_id\nA\n", encoding="utf-8")

        panel = LigandPanel()
        panel.out_dir.set_path(str(out))
        panel._running = True
        panel._prepare_button.setEnabled(False)
        panel._poll()

        assert panel._prepare_button.isEnabled()
        assert not panel._stop_button.isEnabled()


class TestSetupPanelCommand:
    @pytest.fixture
    def panel(self, qapp, tmp_path):
        panel = SetupPanel()
        receptor = tmp_path / "rec.pdbqt"
        receptor.write_text("ATOM\n", encoding="utf-8")
        panel.receptor.set_path(str(receptor))
        panel.ligands.set_path(str(tmp_path / "ligs"))
        panel.run_dir.set_path(str(tmp_path / "run"))
        return panel

    def test_requires_a_box_definition(self, panel):
        assert panel._build_command() is None
        assert "box" in panel._status.text().lower()

    def test_residue_box(self, panel):
        panel.residues.setText("401,405,411")
        command = panel._build_command()

        assert _flag_value(command, "--residues") == "401,405,411"
        assert "--center" not in command

    def test_explicit_box_takes_precedence_over_residues(self, panel):
        panel.residues.setText("401,405")
        panel.center.setText("1,2,3")
        panel.size.setText("20,20,20")
        command = panel._build_command()

        assert _flag_value(command, "--center") == "1,2,3"
        assert "--residues" not in command

    def test_chain_is_included_only_when_given(self, panel):
        panel.residues.setText("401")
        assert "--chain" not in panel._build_command()

        panel.chain.setText("B")
        assert _flag_value(panel._build_command(), "--chain") == "B"

    def test_ad4_engine_requires_maps(self, panel):
        panel.residues.setText("401")
        panel.engine.setCurrentText("ad4")

        assert panel._build_command() is None
        assert "maps" in panel._status.text().lower()

    def test_ad4_engine_with_real_maps_builds(self, panel, tmp_path):
        maps = _fake_maps_dir(tmp_path / "maps")
        panel.residues.setText("401")
        panel.engine.setCurrentText("ad4")
        panel.maps.set_path(str(maps))
        command = panel._build_command()

        assert _flag_value(command, "--engine") == "ad4"
        assert "--maps" in command

    def test_ad4_refuses_a_directory_of_receptors(self, panel, tmp_path):
        """The reported failure: pointing Maps at the receptor directory.

        Left to Vina this surfaces as one identical failure per ligand, so a
        47,000-compound screen fails 47,000 times and explains itself none of
        them. It has to be refused before the run starts.
        """
        receptors = tmp_path / "Receptor"
        receptors.mkdir()
        (receptors / "2OVX.pdbqt").write_text("ATOM\n", encoding="utf-8")
        (receptors / "2OVX_TZ.pdbqt").write_text("ATOM\n", encoding="utf-8")

        panel.residues.setText("401")
        panel.engine.setCurrentText("ad4")
        panel.maps.set_path(str(receptors))

        assert panel._build_command() is None
        message = panel._status.text().lower()
        assert "receptor files" in message
        assert "drydock maps" in message

    def test_maps_are_not_required_for_vina(self, panel):
        panel.residues.setText("401")
        panel.engine.setCurrentText("vina")

        assert panel._build_command() is not None

    def test_missing_inputs_are_named(self, qapp):
        panel = SetupPanel()
        panel.residues.setText("401")
        assert panel._build_command() is None

        message = panel._status.text().lower()
        for expected in ("receptor", "ligand", "run"):
            assert expected in message

    def test_search_settings_are_passed_through(self, panel):
        panel.residues.setText("401")
        panel.exhaustiveness.setValue(16)
        panel.modes.setValue(20)
        panel.seed.setValue(99)
        command = panel._build_command()

        assert _flag_value(command, "--exhaustiveness") == "16"
        assert _flag_value(command, "--modes") == "20"
        assert _flag_value(command, "--seed") == "99"


class TestPathPickerModes:
    """A dialog offering the wrong operation is worse than none.

    Reported: the receptor output field opened a file-*open* dialog, which cannot
    name a file that does not exist yet -- so there was no way to choose an output
    through the interface at all.
    """

    def test_output_fields_use_a_save_dialog(self, qapp):
        from drydock.gui.receptor_panel import ReceptorPanel

        panel = ReceptorPanel()
        assert panel.out._save is True
        assert panel.out._directory is False

    def test_input_fields_use_an_open_dialog(self, qapp):
        from drydock.gui.receptor_panel import ReceptorPanel

        panel = ReceptorPanel()
        assert panel.structure._save is False
        assert panel.structure._directory is False

    def test_directory_fields_stay_directory_pickers(self, qapp):
        from drydock.gui.ligand_panel import LigandPanel

        panel = LigandPanel()
        assert panel.out_dir._directory is True

    def test_basename_field_strips_the_extension_it_will_add(self, qapp):
        """A save dialog naturally produces 'x.pdbqt'; the tool appends it too."""
        from drydock.gui.receptor_panel import ReceptorPanel

        panel = ReceptorPanel()
        panel.out.set_path("/tmp/mmp9.pdbqt")
        assert panel.out.path() == "/tmp/mmp9"

    def test_stripping_is_case_insensitive(self, qapp):
        from drydock.gui.receptor_panel import ReceptorPanel

        panel = ReceptorPanel()
        panel.out.set_path("/tmp/mmp9.PDBQT")
        assert panel.out.path() == "/tmp/mmp9"

    def test_a_plain_basename_is_untouched(self, qapp):
        from drydock.gui.receptor_panel import ReceptorPanel

        panel = ReceptorPanel()
        panel.out.set_path("/tmp/mmp9")
        assert panel.out.path() == "/tmp/mmp9"

    def test_other_extensions_survive(self, qapp):
        """Only the extension the tool adds is removed."""
        from drydock.gui.receptor_panel import ReceptorPanel

        panel = ReceptorPanel()
        panel.out.set_path("/tmp/mmp9.v2")
        assert panel.out.path() == "/tmp/mmp9.v2"


class TestRecprepBasename:
    def test_core_also_strips_a_supplied_extension(self, tmp_path):
        """The CLI can be given 'out.pdbqt' just as easily as the GUI."""
        import inspect as _inspect

        from drydock.core import recprep

        source = _inspect.getsource(recprep.prepare_receptor)
        assert '.with_suffix("")' in source, "basename must be normalised"
