"""Tests for AutoDock4Zn support.

The force-field constants are asserted literally. They are transcribed from a
third-party source, and a typo in one would change every score in a zinc screen
without anything downstream noticing: a mistyped well depth gives plausible
numbers rather than an error. One was transcribed wrong during development,
which is why these are here.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from drydock.core.box import Box
from drydock.core.receptor import CANONICAL_ATOM_TYPES, inspect, normalize_atom_types
from drydock.core.zinc import (
    AD4ZN_PARAMETERS,
    AD4ZN_POTENTIALS,
    GRID_SPACING,
    LIGAND_ATOM_TYPES,
    ZINC_PSEUDO_SCRIPT,
    ZincError,
    add_zinc_pseudo_atoms,
    write_gpf,
)

# The real catalytic zinc site of MMP9 (PDB 2OVX): His401, His405 and His411 in
# full, plus the catalytic zinc.
#
# A hand-built fixture of three isolated nitrogen atoms does not work here, and
# the reason is informative: zinc_pseudo.py infers coordination geometry from
# whole residues, so given bare atoms it finds no coordination and places nothing.
# Testing against real geometry is the only way this path gets exercised.
ZINC_SITE = Path(__file__).parent / "data" / "zinc_site.pdbqt"
ZINC_PDBQT = ZINC_SITE.read_text()

# The catalytic zinc's position in that structure.
CATALYTIC_ZN = (70.219, 18.432, 50.134)


class TestForceFieldConstants:
    """Literal assertions on transcribed parameters."""

    def test_potentials_match_the_published_values(self):
        assert AD4ZN_POTENTIALS == (
            "nbp_r_eps 0.25 3.8581 12 6 NA TZ",
            "nbp_r_eps 2.1  0.6391 12 6 OA Zn",
            "nbp_r_eps 2.25 1.2617 12 6 SA Zn",
            "nbp_r_eps 1.0  0.0    12 6 HD Zn",
            "nbp_r_eps 2.0  0.001  12 6 NA Zn",
            "nbp_r_eps 2.0  0.0493 12 6  N Zn",
        )

    def test_grid_spacing_is_autogrids_default(self):
        """The AD4 function was parameterised at 0.375 A; changing it invalidates it."""
        assert GRID_SPACING == 0.375

    def test_parameter_file_defines_the_pseudo_atom(self):
        """Without a TZ atom_par line, AutoGrid cannot type the pseudo-atoms."""
        assert AD4ZN_PARAMETERS.exists()
        assert "atom_par TZ" in AD4ZN_PARAMETERS.read_text()

    def test_parameter_file_keeps_full_element_coverage(self):
        """TZ was added to the comprehensive set, not swapped for a narrower one."""
        text = AD4ZN_PARAMETERS.read_text()
        assert sum(1 for ln in text.splitlines() if ln.startswith("atom_par")) > 100
        for element in ("atom_par Zn", "atom_par Fe", "atom_par Mg", "atom_par Br"):
            assert element in text

    def test_vendored_script_is_present(self):
        assert ZINC_PSEUDO_SCRIPT.exists()


class TestAtomTypeNormalisation:
    def test_repairs_uppercased_metals(self):
        """zinc_pseudo.py upper-cases metals; Vina rejects the result."""
        line = "HETATM    5 ZN    ZN N 444      70.219  18.432  50.134  1.00 23.42     0.000 ZN\n"
        fixed, changed = normalize_atom_types(line)

        assert changed == 1
        assert fixed[77:79].strip() == "Zn"

    def test_leaves_correct_types_alone(self):
        fixed, changed = normalize_atom_types(ZINC_PDBQT)
        assert changed == 0
        assert fixed == ZINC_PDBQT

    def test_preserves_every_other_column(self):
        line = "HETATM    5 ZN    ZN N 444      70.219  18.432  50.134  1.00 23.42     0.000 ZN\n"
        fixed, _ = normalize_atom_types(line)

        assert fixed[:77] == line[:77], "coordinates and charges must not move"
        assert len(fixed) == len(line)

    @pytest.mark.parametrize("element", ["ZN", "MG", "FE", "MN", "CL", "BR"])
    def test_known_metals_and_halogens_are_canonicalised(self, element):
        assert CANONICAL_ATOM_TYPES[element] == element.capitalize()

    def test_single_letter_types_are_unchanged(self):
        for t in ("C", "N", "P", "S", "F", "I", "A"):
            assert CANONICAL_ATOM_TYPES[t] == t

    def test_ignores_non_atom_lines(self):
        text = "REMARK this mentions ZN but is not an atom record\n"
        fixed, changed = normalize_atom_types(text)
        assert changed == 0
        assert fixed == text


class TestAddZincPseudoAtoms:
    def test_places_a_pseudo_atom_at_the_vacant_site(self, tmp_path):
        receptor = tmp_path / "zn.pdbqt"
        receptor.write_text(ZINC_PDBQT, encoding="utf-8")

        result = add_zinc_pseudo_atoms(receptor, tmp_path / "zn_TZ.pdbqt")

        assert result.n_zinc == 1
        assert result.n_pseudo_atoms >= 1
        assert result.receptor_tz.exists()

    def test_pseudo_atoms_sit_at_the_coordination_distance(self, tmp_path):
        """AutoDock4Zn places TZ 2.0 A from the zinc, which is the whole point."""
        import math

        receptor = tmp_path / "zn.pdbqt"
        receptor.write_text(ZINC_PDBQT, encoding="utf-8")
        result = add_zinc_pseudo_atoms(receptor, tmp_path / "zn_TZ.pdbqt")

        zinc = CATALYTIC_ZN
        distances = [
            math.dist(zinc, (float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
            for ln in result.receptor_tz.read_text().splitlines()
            if ln.startswith(("ATOM", "HETATM")) and ln[77:79].strip() == "TZ"
        ]
        assert distances
        for d in distances:
            assert d == pytest.approx(2.0, abs=0.05)

    def test_output_types_survive_as_valid_autodock_types(self, tmp_path):
        """The output must be usable with Vina, not only with AutoGrid."""
        receptor = tmp_path / "zn.pdbqt"
        receptor.write_text(ZINC_PDBQT, encoding="utf-8")
        result = add_zinc_pseudo_atoms(receptor, tmp_path / "zn_TZ.pdbqt")

        report = inspect(result.receptor_tz)
        assert "Zn" in report.atom_types, "metal casing must be repaired"
        assert "ZN" not in report.atom_types
        assert report.has_zinc_pseudo_atoms

    def test_receptor_without_zinc_is_refused_with_a_reason(self, tmp_path):
        receptor = tmp_path / "nozn.pdbqt"
        receptor.write_text(
            "ATOM      1  CA  ALA B 110      10.000  10.000  10.000  1.00 20.00     0.000 C\n",
            encoding="utf-8",
        )
        with pytest.raises(ZincError, match="no zinc"):
            add_zinc_pseudo_atoms(receptor, tmp_path / "out.pdbqt")

    def test_missing_receptor_is_reported(self, tmp_path):
        with pytest.raises(ZincError, match="not found"):
            add_zinc_pseudo_atoms(tmp_path / "absent.pdbqt")


class TestWriteGpf:
    @pytest.fixture
    def tz_receptor(self, tmp_path):
        receptor = tmp_path / "zn.pdbqt"
        receptor.write_text(ZINC_PDBQT, encoding="utf-8")
        return add_zinc_pseudo_atoms(receptor, tmp_path / "zn_TZ.pdbqt").receptor_tz

    def test_includes_every_ad4zn_potential(self, tz_receptor, tmp_path):
        box = Box((70.0, 18.0, 50.0), (20.0, 20.0, 20.0))
        gpf = write_gpf(tz_receptor, box, tmp_path / "g.gpf")
        text = gpf.read_text()

        for potential in AD4ZN_POTENTIALS:
            assert potential in text

    def test_grid_points_are_odd_so_the_box_is_centred(self, tz_receptor, tmp_path):
        """An even count would centre the grid between two points, not on one."""
        box = Box((70.0, 18.0, 50.0), (20.0, 21.0, 22.0))
        gpf = write_gpf(tz_receptor, box, tmp_path / "g.gpf")

        npts_line = next(ln for ln in gpf.read_text().splitlines() if ln.startswith("npts"))
        counts = [int(v) for v in npts_line.split()[1:]]
        assert all(n % 2 == 0 for n in counts), "AutoGrid npts are point intervals"

    def test_grid_centre_matches_the_box(self, tz_receptor, tmp_path):
        box = Box((1.5, 2.5, 3.5), (20.0, 20.0, 20.0))
        gpf = write_gpf(tz_receptor, box, tmp_path / "g.gpf")

        line = next(ln for ln in gpf.read_text().splitlines() if ln.startswith("gridcenter"))
        assert [float(v) for v in line.split()[1:]] == [1.5, 2.5, 3.5]

    def test_receptor_types_include_the_pseudo_atom(self, tz_receptor, tmp_path):
        box = Box((70.0, 18.0, 50.0), (20.0, 20.0, 20.0))
        gpf = write_gpf(tz_receptor, box, tmp_path / "g.gpf")

        line = next(ln for ln in gpf.read_text().splitlines() if ln.startswith("receptor_types"))
        assert "TZ" in line.split()
        assert "Zn" in line.split()

    def test_ligand_types_cover_what_meeko_emits(self, tz_receptor, tmp_path):
        """A missing map is a hard failure at docking time, not a skipped atom."""
        box = Box((70.0, 18.0, 50.0), (20.0, 20.0, 20.0))
        gpf = write_gpf(tz_receptor, box, tmp_path / "g.gpf")

        line = next(ln for ln in gpf.read_text().splitlines() if ln.startswith("ligand_types"))
        declared = set(line.split()[1:])
        assert declared == set(LIGAND_ATOM_TYPES)
        for essential in ("C", "A", "N", "NA", "OA", "HD", "SA"):
            assert essential in declared

    def test_requests_one_map_per_ligand_type(self, tz_receptor, tmp_path):
        box = Box((70.0, 18.0, 50.0), (20.0, 20.0, 20.0))
        gpf = write_gpf(tz_receptor, box, tmp_path / "g.gpf")

        maps = [ln for ln in gpf.read_text().splitlines() if ln.startswith("map ")]
        assert len(maps) == len(LIGAND_ATOM_TYPES)

    def test_points_at_the_parameter_file(self, tz_receptor, tmp_path):
        box = Box((70.0, 18.0, 50.0), (20.0, 20.0, 20.0))
        gpf = write_gpf(tz_receptor, box, tmp_path / "g.gpf")
        assert "AD4Zn.dat" in gpf.read_text()

    def test_accepts_a_custom_parameter_file(self, tz_receptor, tmp_path):
        custom = tmp_path / "mine.dat"
        custom.write_text("atom_par TZ 1.0 0.0 0.0 0.0 0.0 0.0 0 -1 -1 0\n", encoding="utf-8")
        gpf = write_gpf(
            tz_receptor,
            Box((70.0, 18.0, 50.0), (20.0, 20.0, 20.0)),
            tmp_path / "g.gpf",
            parameter_file=custom,
        )
        assert "mine.dat" in gpf.read_text()
