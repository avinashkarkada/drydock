"""Tests for receptor preparation.

These do not test that Meeko works, which is Meeko's business. They test that
its output is checked before being handed on, and that its failures arrive as
something a user can act on. A receptor that prepares "successfully" with no
hydrogen-bond donors is the failure that prompted all this, and it raises
nothing anywhere.
"""

from __future__ import annotations

import pytest

from drydock.core.recprep import (
    CIF_FORMATS,
    NATIVE_FORMATS,
    ReceptorPrepError,
    _deleted_residues,
    _explain_failure,
    to_pdb,
)

CIF = "/home/aka/Desktop/MS/Dr.Pavan_G/Virtual_screen/2OVX.cif"


class TestFormatHandling:
    def test_native_formats_pass_through_untouched(self, tmp_path):
        pdb = tmp_path / "x.pdb"
        pdb.write_text("ATOM\n", encoding="utf-8")
        assert to_pdb(pdb) == pdb

    def test_unsupported_format_is_refused_with_the_alternatives(self, tmp_path):
        bad = tmp_path / "x.xyz"
        bad.write_text("", encoding="utf-8")

        with pytest.raises(ReceptorPrepError, match="unsupported"):
            to_pdb(bad)

    def test_format_sets_do_not_overlap(self):
        assert not (NATIVE_FORMATS & CIF_FORMATS)

    @pytest.mark.skipif(not __import__("pathlib").Path(CIF).exists(), reason="needs 2OVX.cif")
    def test_converts_mmcif_to_pdb(self, tmp_path):
        out = to_pdb(CIF, tmp_path / "converted.pdb")

        assert out.exists()
        text = out.read_text()
        assert text.startswith(("HEADER", "CRYST", "ATOM", "REMARK", "MODEL", "EXPDTA"))
        assert sum(1 for ln in text.splitlines() if ln.startswith(("ATOM", "HETATM"))) > 100

    @pytest.mark.skipif(not __import__("pathlib").Path(CIF).exists(), reason="needs 2OVX.cif")
    def test_conversion_keeps_the_metals(self, tmp_path):
        """Losing a zinc during format conversion would be silent and fatal."""
        out = to_pdb(CIF, tmp_path / "converted.pdb")
        zincs = [
            ln
            for ln in out.read_text().splitlines()
            if ln.startswith("HETATM") and "ZN" in ln[17:20].upper()
        ]
        assert len(zincs) == 2

    def test_missing_file_is_reported(self, tmp_path):
        from drydock.core.recprep import prepare_receptor

        with pytest.raises(ReceptorPrepError, match="not found"):
            prepare_receptor(tmp_path / "absent.pdb", tmp_path / "out")


class TestFailureExplanations:
    """Meeko's failures are correct but arrive as tracebacks."""

    def test_missing_atoms_names_the_flag_that_helps(self):
        message = _explain_failure("RuntimeError: residue has missing atoms")
        assert "--allow-bad-residues" in message

    def test_altloc_failure_names_the_flag(self):
        message = _explain_failure("ValueError: ambiguous altloc for residue 42")
        assert "--altloc" in message

    def test_template_failure_suggests_what_it_usually_means(self):
        message = _explain_failure("no template found for residue LIG")
        assert "modified residue" in message or "template" in message
        assert "--delete-residues" in message

    def test_unrecognised_failure_still_surfaces_the_output(self):
        message = _explain_failure("Segmentation fault somewhere deep")
        assert "Segmentation fault" in message

    def test_explanations_are_bounded(self):
        assert len(_explain_failure("x" * 10_000)) < 1200


class TestDeletedResidueReporting:
    def test_extracts_deleted_residues(self):
        """A silently deleted active-site residue changes the pocket."""
        output = "some noise\nDeleted residue A:350 (missing atoms)\nmore noise\n"
        assert any("A:350" in line for line in _deleted_residues(output))

    def test_reports_nothing_when_nothing_was_deleted(self):
        assert _deleted_residues("prepared successfully\n") == ()


@pytest.mark.skipif(not __import__("pathlib").Path(CIF).exists(), reason="needs 2OVX.cif")
class TestEndToEnd:
    """Against the real structure, since that is where the failure came from."""

    @pytest.fixture(scope="class")
    def prepared(self, tmp_path_factory):
        from drydock.core.recprep import prepare_receptor

        out = tmp_path_factory.mktemp("recprep") / "receptor"
        return prepare_receptor(CIF, out)

    def test_produces_a_pdbqt(self, prepared):
        assert prepared.receptor_pdbqt.exists()
        assert prepared.receptor_pdbqt.suffix == ".pdbqt"

    def test_records_that_it_converted_from_cif(self, prepared):
        assert prepared.converted_from is not None
        assert prepared.converted_from.suffix == ".cif"

    def test_adds_polar_hydrogens(self, prepared):
        """The whole point. The input receptor had none."""
        assert prepared.report.n_polar_hydrogens > 100

    def test_uses_canonical_metal_casing(self, prepared):
        """'ZN' would be rejected by Vina at load time."""
        assert "Zn" in prepared.report.atom_types
        assert "ZN" not in prepared.report.atom_types

    def test_keeps_both_zincs(self, prepared):
        zinc_atoms = prepared.report.atom_types.get("Zn", 0)
        assert zinc_atoms == 2

    def test_result_passes_its_own_check(self, prepared):
        assert prepared.ok, f"unexpected problems: {prepared.report.problems}"

    def test_intermediate_pdb_is_cleaned_up(self, prepared):
        leftovers = list(prepared.receptor_pdbqt.parent.glob("*.converted.pdb"))
        assert leftovers == []

    def test_output_loads_in_vina(self, prepared):
        """The end-to-end assertion: the engine accepts what preparation produced."""
        from vina import Vina

        v = Vina(sf_name="vina", cpu=1, seed=1, verbosity=0)
        v.set_receptor(str(prepared.receptor_pdbqt))  # must not raise
