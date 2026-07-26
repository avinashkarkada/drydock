"""Tests for ligand preparation.

Includes regressions for two bugs that only appeared against real data:

* Torsion counts read from the wrong PDBQT record, reporting 0 rotatable bonds
  for every compound in the library.
* Scrubber mutating the molecule it is handed -- stripping its conformer even
  when it fails -- so the geometry fallback was given an already-emptied
  molecule and reported "no 3D conformer produced" for compounds whose input
  coordinates were fine.
"""

from __future__ import annotations

import pytest

from drydock.core.descriptors import compute, ligand_efficiency
from drydock.core.library import Record
from drydock.core.ligprep import (
    PreparedLigand,
    PrepConfig,
    PrepFailure,
    _count_torsions,
    _readable_error,
    prepare_one,
)
from drydock.core.seeds import ligand_seed


def _smiles_record(ligand_id: str, smiles: str) -> Record:
    return Record(ligand_id, ligand_id, smiles, "smi", 0)


class TestDescriptors:
    def test_computes_expected_values_for_ethanol(self):
        from rdkit import Chem

        desc = compute(Chem.MolFromSmiles("CCO"))
        assert desc.formula == "C2H6O"
        assert desc.heavy_atoms == 3
        assert desc.mw == pytest.approx(46.07, abs=0.01)
        assert desc.hbd == 1

    def test_hydrogens_do_not_change_heavy_atom_count(self):
        """Descriptors must describe the compound, not the file's H treatment."""
        from rdkit import Chem

        bare = Chem.MolFromSmiles("CCO")
        with_h = Chem.AddHs(bare)
        assert compute(bare).heavy_atoms == compute(with_h).heavy_atoms == 3

    def test_formal_charge_is_reported(self):
        from rdkit import Chem

        assert compute(Chem.MolFromSmiles("CC(=O)[O-]")).formal_charge == -1


class TestLigandEfficiency:
    def test_divides_affinity_by_heavy_atoms(self):
        assert ligand_efficiency(-10.0, 25) == pytest.approx(-0.4)

    @pytest.mark.parametrize(("affinity", "heavy"), [(None, 25), (-10.0, 0), (-10.0, None)])
    def test_missing_inputs_give_none_not_an_error(self, affinity, heavy):
        assert ligand_efficiency(affinity, heavy) is None

    def test_penalises_size(self):
        """The whole point: a bigger molecule needs a better score to rank."""
        small = ligand_efficiency(-8.0, 20)
        large = ligand_efficiency(-9.0, 40)
        assert small < large, "the smaller compound should be the more efficient"


class TestTorsionCounting:
    def test_reads_the_torsdof_record(self):
        """Meeko writes TORSDOF at the end of the file, not a header REMARK."""
        pdbqt = "ROOT\nATOM      1  C\nENDROOT\nBRANCH   1   2\nENDBRANCH   1   2\nTORSDOF 7\n"
        assert _count_torsions(pdbqt) == 7

    def test_falls_back_to_counting_branches(self):
        pdbqt = "ROOT\nENDROOT\nBRANCH   1   2\nBRANCH   2   3\nENDBRANCH   2   3\n"
        assert _count_torsions(pdbqt) == 2

    def test_accepts_the_mgltools_header_format(self):
        """PDBQTs from prepare_ligand4.py are still widely in circulation."""
        pdbqt = "REMARK  18 active torsions:\nROOT\nATOM      1  C\n"
        assert _count_torsions(pdbqt) == 18

    def test_rigid_molecule_has_no_torsions(self):
        assert _count_torsions("ROOT\nATOM      1  C\nENDROOT\nTORSDOF 0\n") == 0


class TestReadableError:
    def test_translates_rdkit_embedding_counters(self):
        """RDKit reports embedding failure as an internal stage tally."""
        raw = Exception("{'INITIAL_COORDS': 230, 'FIRST_MINIMIZATION': 1}")
        message = _readable_error(raw)
        assert "3D embedding failed" in message

    def test_passes_ordinary_messages_through(self):
        assert _readable_error(ValueError("something specific")) == "something specific"

    def test_truncates_runaway_messages(self):
        assert len(_readable_error(ValueError("x" * 5000))) <= 200


class TestPrepareOne:
    def test_prepares_a_simple_molecule(self, tmp_path):
        result = prepare_one(_smiles_record("ethanol", "CCO"), PrepConfig(seed=1), tmp_path)

        assert isinstance(result, PreparedLigand)
        assert result.pdbqt_files == ("ethanol.pdbqt",)
        assert (tmp_path / "ethanol.pdbqt").exists()
        assert result.descriptors.formula == "C2H6O"

    def test_counts_torsions_of_a_flexible_molecule(self, tmp_path):
        """The regression: this reported 0 for every compound in the library."""
        result = prepare_one(
            _smiles_record("chain", "CCCCCCCCCC"), PrepConfig(seed=1), tmp_path
        )

        assert isinstance(result, PreparedLigand)
        assert result.torsions > 0, "a decane chain has rotatable bonds"

    def test_unparseable_input_is_a_failure_not_an_exception(self, tmp_path):
        result = prepare_one(
            _smiles_record("junk", "not a smiles!!"), PrepConfig(seed=1), tmp_path
        )

        assert isinstance(result, PrepFailure)
        assert result.stage == "parse"
        assert result.ligand_id == "junk"

    def test_failure_keeps_both_identifiers(self, tmp_path):
        """A failure must be reportable against a name the user recognises."""
        record = Record("CPD_7", "CPD", "garbage((", "smi", 3)
        result = prepare_one(record, PrepConfig(seed=1), tmp_path)

        assert isinstance(result, PrepFailure)
        assert (result.ligand_id, result.compound_id) == ("CPD_7", "CPD")

    def test_same_seed_gives_byte_identical_output(self, tmp_path):
        a, b = tmp_path / "a", tmp_path / "b"
        record = _smiles_record("mol", "CC(=O)Oc1ccccc1C(=O)O")

        prepare_one(record, PrepConfig(seed=7), a)
        prepare_one(record, PrepConfig(seed=7), b)

        assert (a / "mol.pdbqt").read_text() == (b / "mol.pdbqt").read_text()

    def test_max_torsions_filter_rejects_floppy_molecules(self, tmp_path):
        record = _smiles_record("floppy", "CCCCCCCCCCCCCCCCCCCC")
        result = prepare_one(record, PrepConfig(seed=1, max_torsions=3), tmp_path)

        assert isinstance(result, PrepFailure)
        assert result.stage == "filter"
        assert "max_torsions" in result.error

    def test_manifest_row_is_flat(self, tmp_path):
        result = prepare_one(_smiles_record("ethanol", "CCO"), PrepConfig(seed=1), tmp_path)
        row = result.to_row()

        assert row["ligand_id"] == "ethanol"
        assert row["compound_id"] == "ethanol"
        assert isinstance(row["mw"], float)
        assert "smiles" in row


class TestGeometryFallback:
    """Regression tests for Scrubber mutating the molecule it is given."""

    def test_embedding_failure_falls_back_to_input_geometry(self, tmp_path, monkeypatch):
        from rdkit import Chem

        from drydock.core import ligprep

        # A molecule that already carries 3D coordinates, as SDF input does.
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        from rdkit.Chem import AllChem

        AllChem.EmbedMolecule(mol, randomSeed=1)
        block = Chem.MolToMolBlock(mol)
        record = Record("has3d", "has3d", block, "sdf", 0)

        tools = ligprep._get_tools(PrepConfig(seed=1))
        real_scrub = tools.scrub

        def always_fails(_mol):
            # Reproduce the real failure precisely: Scrubber strips the caller's
            # conformer on the way out, even when it raises.
            _mol.RemoveAllConformers()
            raise RuntimeError("{'INITIAL_COORDS': 230}")

        monkeypatch.setattr(tools, "scrub", always_fails)
        try:
            result = prepare_one(record, PrepConfig(seed=1), tmp_path)
        finally:
            tools.scrub = real_scrub

        assert isinstance(result, PreparedLigand), "input geometry should have rescued it"
        assert any("input geometry" in w for w in result.warnings)

    def test_smiles_input_cannot_fall_back(self, tmp_path, monkeypatch):
        """There is no input geometry to fall back to, so failing is correct."""
        from drydock.core import ligprep

        tools = ligprep._get_tools(PrepConfig(seed=1))
        real_scrub = tools.scrub

        def always_fails(_mol):
            raise RuntimeError("{'INITIAL_COORDS': 230}")

        monkeypatch.setattr(tools, "scrub", always_fails)
        try:
            result = prepare_one(_smiles_record("flat", "CCO"), PrepConfig(seed=1), tmp_path)
        finally:
            tools.scrub = real_scrub

        assert isinstance(result, PrepFailure)
        assert result.stage == "protonate"


class TestSeeds:
    def test_is_stable_across_calls(self):
        assert ligand_seed(42, "CMNPD1") == ligand_seed(42, "CMNPD1")

    def test_differs_by_ligand_and_by_global_seed(self):
        assert ligand_seed(42, "A") != ligand_seed(42, "B")
        assert ligand_seed(1, "A") != ligand_seed(2, "A")

    def test_does_not_depend_on_scheduling_order(self):
        """Derived from the identifier, so parallelism cannot perturb it."""
        forwards = [ligand_seed(9, f"L{i}") for i in range(50)]
        backwards = [ligand_seed(9, f"L{i}") for i in reversed(range(50))]
        assert forwards == list(reversed(backwards))

    def test_fits_in_a_signed_32_bit_range(self):
        """Vina rejects anything wider."""
        for i in range(500):
            assert 0 <= ligand_seed(i, f"lig{i}") <= 2**31 - 1
