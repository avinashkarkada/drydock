"""Tests for the search box and receptor checking.

The receptor tests encode two failure modes found in a real prepared receptor,
both of which produce a valid file that then ruins a screen:

* No polar hydrogens, so the protein has no hydrogen-bond donors at all.
* Atom types in the wrong case (``ZN`` for ``Zn``), which Vina rejects with a
  C++ overload error that reads like a bug in the caller.
"""

from __future__ import annotations

import pytest

from drydock.core.box import LARGE_BOX_VOLUME, Box
from drydock.core.receptor import box_from_residues, inspect, read_pdbqt, select_residues

# Two residues, four atoms, spanning a known volume.
PDBQT = """\
REMARK   4 XXXX COMPLIES WITH FORMAT V. 2.0
ATOM      1  N   PHE B 110      10.000  10.000  10.000  1.00 58.87    -0.229 N
ATOM      2  CA  PHE B 110      12.000  10.000  10.000  1.00 60.17     0.186 C
ATOM      3  CB  PHE B 110      12.000  14.000  10.000  1.00 59.73     0.034 C
ATOM      4  N   LYS B 111      10.000  10.000  16.000  1.00 32.24    -0.229 N
ATOM      5  HN  LYS B 111      10.500  10.500  16.500  1.00  0.00     0.275 HD
HETATM    6 ZN    ZN N 444      11.000  11.000  12.000  1.00 23.42     0.000 Zn
"""


@pytest.fixture
def receptor(tmp_path):
    path = tmp_path / "rec.pdbqt"
    path.write_text(PDBQT, encoding="utf-8")
    return path


class TestBox:
    def test_geometry(self):
        box = Box((0.0, 0.0, 0.0), (10.0, 20.0, 30.0))
        assert box.volume == 6000
        assert box.minimum == (-5.0, -10.0, -15.0)
        assert box.maximum == (5.0, 10.0, 15.0)

    def test_rejects_non_positive_dimensions(self):
        with pytest.raises(ValueError, match="positive"):
            Box((0.0, 0.0, 0.0), (10.0, 0.0, 10.0))

    def test_contains(self):
        box = Box((0.0, 0.0, 0.0), (10.0, 10.0, 10.0))
        assert box.contains((0.0, 0.0, 0.0))
        assert box.contains((5.0, 5.0, 5.0)), "boundary counts as inside"
        assert not box.contains((6.0, 0.0, 0.0))

    def test_from_atoms_encloses_everything_with_padding(self):
        coords = [(0.0, 0.0, 0.0), (10.0, 0.0, 0.0), (0.0, 10.0, 0.0)]
        box = Box.from_atoms(coords, padding=2.0)

        assert box.size[0] == pytest.approx(14.0)
        for point in coords:
            assert box.contains(point)

    def test_from_atoms_cubic_uses_the_largest_dimension(self):
        box = Box.from_atoms([(0.0, 0.0, 0.0), (10.0, 2.0, 2.0)], padding=0.0, cubic=True)
        assert box.size == (10.0, 10.0, 10.0)

    def test_from_atoms_rejects_empty_input(self):
        with pytest.raises(ValueError, match="no atoms"):
            Box.from_atoms([])

    def test_warns_about_large_boxes(self):
        """Vina's own threshold: search quality degrades in large volumes."""
        box = Box((0.0, 0.0, 0.0), (40.0, 40.0, 40.0))
        assert box.volume > LARGE_BOX_VOLUME
        assert any("volume" in w for w in box.warnings())

    def test_warns_about_boxes_too_small_for_a_ligand(self):
        assert any("small" in w for w in Box((0.0, 0.0, 0.0), (8.0, 20.0, 20.0)).warnings())

    def test_reasonable_box_has_no_warnings(self):
        assert Box((0.0, 0.0, 0.0), (22.0, 22.0, 22.0)).warnings() == []

    def test_vina_config_roundtrip(self):
        box = Box((1.5, 2.5, 3.5), (20.0, 21.0, 22.0))
        text = box.to_vina_config()

        parsed = {}
        for line in text.splitlines():
            if "=" in line:
                key, value = line.split("=")
                parsed[key.strip()] = float(value)
        assert Box.from_config(parsed) == box

    def test_from_config_accepts_the_list_form(self):
        box = Box.from_config({"center": [1, 2, 3], "size": [10, 11, 12]})
        assert box.center == (1.0, 2.0, 3.0)

    def test_from_config_reports_what_is_missing(self):
        with pytest.raises(ValueError, match="missing"):
            Box.from_config({"center_x": 1.0})


class TestReceptorReading:
    def test_reads_atoms_by_column(self, receptor):
        atoms = read_pdbqt(receptor)
        assert len(atoms) == 6
        assert atoms[0].residue == "PHE"
        assert atoms[0].coordinates == (10.0, 10.0, 10.0)
        assert atoms[-1].atom_type == "Zn"
        assert atoms[-1].is_metal

    def test_column_parsing_survives_wide_coordinates(self, tmp_path):
        """Whitespace splitting breaks once adjacent fields run together.

        Built to exact PDB column positions so the coordinates abut with no
        separating space, which is legal and happens in real structures as soon
        as a coordinate needs four digits or a B-factor reaches 100.
        """
        # Each coordinate fills its 8-column field exactly, so the three run
        # together with no separating whitespace at all.
        line = (
            f"{'ATOM':<6}{1:>5} {' N  ':<4}{'':1}{'PHE':>3} {'B':1}{110:>4}    "
            f"{-100.123:>8.3f}{-100.456:>8.3f}{-999.789:>8.3f}"
            f"{1.00:>6.2f}{100.00:>6.2f}    {-0.229:>6.3f} {'N':<2}\n"
        )
        assert " " not in line[30:54], "coordinates must abut for this to be a real test"

        path = tmp_path / "wide.pdbqt"
        path.write_text(line, encoding="utf-8")

        atoms = read_pdbqt(path)
        assert len(atoms) == 1
        assert atoms[0].coordinates == (-100.123, -100.456, -999.789)


class TestReceptorInspection:
    def test_reports_composition(self, receptor):
        report = inspect(receptor)
        assert report.n_atoms == 6
        assert report.metals == ("ZN",)
        assert report.residue_range == (110, 444)
        assert report.n_polar_hydrogens == 1

    def test_missing_polar_hydrogens_is_a_problem(self, tmp_path):
        """The real failure: a receptor with no hydrogen-bond donors at all."""
        text = "\n".join(ln for ln in PDBQT.splitlines() if not ln.endswith("HD"))
        path = tmp_path / "noh.pdbqt"
        path.write_text(text + "\n", encoding="utf-8")

        report = inspect(path)
        assert not report.ok
        assert any("polar hydrogen" in p for p in report.problems)

    def test_zinc_without_pseudo_atoms_is_a_note_not_a_problem(self, receptor):
        report = inspect(receptor)
        assert report.ok
        assert any("TZ pseudo-atoms" in n for n in report.notes)

    def test_empty_file_is_reported_clearly(self, tmp_path):
        path = tmp_path / "empty.pdbqt"
        path.write_text("", encoding="utf-8")
        assert any("no ATOM" in p for p in inspect(path).problems)

    def test_pdb_without_types_is_flagged(self, tmp_path):
        path = tmp_path / "plain.pdb"
        path.write_text(
            "ATOM      1  N   PHE B 110      10.000  10.000  10.000  1.00 58.87\n",
            encoding="utf-8",
        )
        report = inspect(path)
        assert any("PDBQT" in p for p in report.problems)


class TestResidueSelection:
    def test_selects_by_residue_number(self, receptor):
        atoms = read_pdbqt(receptor)
        assert len(select_residues(atoms, [110])) == 3

    def test_sidechains_only_drops_backbone(self, receptor):
        atoms = read_pdbqt(receptor)
        selected = select_residues(atoms, [110], sidechains_only=True)
        assert {a.name for a in selected} == {"CB"}

    def test_chain_filter(self, receptor):
        atoms = read_pdbqt(receptor)
        assert select_residues(atoms, [110], chain="Z") == []


class TestBoxFromResidues:
    def test_encloses_the_named_residues(self, receptor):
        box, atoms = box_from_residues(receptor, [110, 111], padding=2.0)
        assert len(atoms) == 5
        for atom in atoms:
            assert box.contains(atom.coordinates)

    def test_missing_residue_is_an_error_not_a_smaller_box(self, receptor):
        """A mistyped residue number would otherwise shrink the box silently."""
        with pytest.raises(ValueError, match="not present"):
            box_from_residues(receptor, [110, 999], padding=2.0)

    def test_error_mentions_numbering(self, receptor):
        with pytest.raises(ValueError, match="numbering"):
            box_from_residues(receptor, [999], padding=2.0)

    def test_no_matching_residues_at_all(self, receptor):
        with pytest.raises(ValueError, match="none of residues"):
            box_from_residues(receptor, [500, 501], padding=2.0)

    def test_padding_widens_the_box(self, receptor):
        tight, _ = box_from_residues(receptor, [110, 111], padding=1.0)
        loose, _ = box_from_residues(receptor, [110, 111], padding=6.0)
        assert loose.volume > tight.volume
        assert loose.center == tight.center, "padding must not move the centre"
