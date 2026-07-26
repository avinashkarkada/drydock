"""Tests for streaming library readers.

Heavily weighted towards identifier handling, because that is where a library
reader silently destroys data. The motivating case is real: CMNPD 1.0 has 47,451
records under 25,224 distinct ``COMPOUND_ID`` values, and a reader that treats
the identifier as a filename keeps one stereoisomer of each and loses the rest
without reporting anything.
"""

from __future__ import annotations

import gzip

import pytest

from drydock.core.library import (
    LibraryFormatError,
    count_records,
    detect_format,
    iter_library,
    sanitize_id,
    survey,
)

# A minimal but valid SDF: two records, ID in a COMPOUND_ID property.
SDF_TEXT = """\

  test

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0
M  END
> <COMPOUND_ID>
CPD1

$$$$

  test

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 O   0  0
M  END
> <COMPOUND_ID>
CPD2

$$$$
"""


def _write(tmp_path, name: str, text: str):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


class TestFormatDetection:
    @pytest.mark.parametrize(
        ("name", "expected"),
        [
            ("lib.sdf", "sdf"),
            ("lib.SDF", "sdf"),
            ("lib.sdf.gz", "sdf"),
            ("lib.smi", "smi"),
            ("lib.smiles", "smi"),
            ("lib.mol2", "mol2"),
            ("lib.mol2.gz", "mol2"),
        ],
    )
    def test_detects_by_extension(self, name, expected):
        assert detect_format(name) == expected

    def test_unknown_extension_is_an_error(self):
        with pytest.raises(LibraryFormatError):
            detect_format("compounds.xyz")


class TestSanitizeId:
    @pytest.mark.parametrize(
        ("raw", "expected"),
        [
            ("CMNPD1", "CMNPD1"),
            ("  CMNPD1  ", "CMNPD1"),
            ("a/b", "a_b"),
            ("../escape", ".._escape"),
            ("name with spaces", "name_with_spaces"),
            ("﻿CMNPD1", "CMNPD1"),
        ],
    )
    def test_makes_ids_filename_safe(self, raw, expected):
        assert sanitize_id(raw) == expected

    def test_path_separators_cannot_survive(self):
        """An ID that kept a slash would write outside the output directory."""
        assert "/" not in sanitize_id("evil/../../etc/passwd")


class TestSdfReading:
    def test_reads_records_and_ids(self, tmp_path):
        records = list(iter_library(_write(tmp_path, "l.sdf", SDF_TEXT)))
        assert [r.ligand_id for r in records] == ["CPD1", "CPD2"]
        assert [r.index for r in records] == [0, 1]
        assert all(r.fmt == "sdf" for r in records)

    def test_reads_gzipped_libraries(self, tmp_path):
        path = tmp_path / "l.sdf.gz"
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            fh.write(SDF_TEXT)

        assert [r.ligand_id for r in iter_library(path)] == ["CPD1", "CPD2"]

    def test_final_record_without_terminator_is_kept(self, tmp_path):
        truncated = SDF_TEXT.rstrip()[: -len("$$$$")]
        records = list(iter_library(_write(tmp_path, "l.sdf", truncated)))
        assert [r.ligand_id for r in records] == ["CPD1", "CPD2"]

    def test_byte_order_mark_is_not_an_identifier(self, tmp_path):
        """A BOM leads the real CMNPD file and was read as the first ID."""
        records = list(iter_library(_write(tmp_path, "l.sdf", "﻿" + SDF_TEXT)))
        assert records[0].ligand_id == "CPD1"

    def test_explicit_id_field_wins(self, tmp_path):
        text = SDF_TEXT.replace("> <COMPOUND_ID>", "> <CATALOG_ID>")
        records = list(iter_library(_write(tmp_path, "l.sdf", text), id_field="CATALOG_ID"))
        assert [r.ligand_id for r in records] == ["CPD1", "CPD2"]

    def test_falls_back_to_a_positional_id(self, tmp_path):
        text = SDF_TEXT.replace("> <COMPOUND_ID>\nCPD1\n", "").replace(
            "> <COMPOUND_ID>\nCPD2\n", ""
        )
        records = list(iter_library(_write(tmp_path, "l.sdf", text)))
        assert all(r.ligand_id.startswith("mol") for r in records)

    def test_title_line_is_used_when_no_property_matches(self, tmp_path):
        # Line 0 of an SDF record is the title; line 1 is the program banner.
        text = (
            "MYNAME\n  test\n\n"
            "  1  0  0  0  0  0            999 V2000\n"
            "    0.0000    0.0000    0.0000 C   0  0\nM  END\n$$$$\n"
        )
        assert list(iter_library(_write(tmp_path, "l.sdf", text)))[0].ligand_id == "MYNAME"

    def test_property_id_beats_the_title_line(self, tmp_path):
        """Bulk-generated SDFs put banners in the title; a tag is deliberate."""
        text = (
            "SciTegic3D\n  test\n\n"
            "  1  0  0  0  0  0            999 V2000\n"
            "    0.0000    0.0000    0.0000 C   0  0\nM  END\n"
            "> <COMPOUND_ID>\nREAL_ID\n\n$$$$\n"
        )
        assert list(iter_library(_write(tmp_path, "l.sdf", text)))[0].ligand_id == "REAL_ID"


class TestDuplicateIdentifiers:
    """The stereoisomer case: identifiers repeat and every record must survive."""

    def _library(self, tmp_path, ids: list[str]):
        blocks = []
        for cid in ids:
            blocks.append(
                "\n  test\n\n  1  0  0  0  0  0            999 V2000\n"
                "    0.0000    0.0000    0.0000 C   0  0\nM  END\n"
                f"> <COMPOUND_ID>\n{cid}\n\n$$$$\n"
            )
        return _write(tmp_path, "l.sdf", "".join(blocks))

    def test_repeats_get_distinct_ligand_ids(self, tmp_path):
        path = self._library(tmp_path, ["A", "A", "A", "B"])
        records = list(iter_library(path))

        assert [r.ligand_id for r in records] == ["A", "A_2", "A_3", "B"]
        assert [r.compound_id for r in records] == ["A", "A", "A", "B"]

    def test_no_record_is_lost(self, tmp_path):
        path = self._library(tmp_path, ["X"] * 64)
        records = list(iter_library(path))

        assert len(records) == 64
        assert len({r.ligand_id for r in records}) == 64, "ids must be unique"
        assert {r.compound_id for r in records} == {"X"}, "grouping must be preserved"

    def test_is_duplicate_id_flags_renamed_records(self, tmp_path):
        records = list(iter_library(self._library(tmp_path, ["A", "A"])))
        assert records[0].is_duplicate_id is False
        assert records[1].is_duplicate_id is True

    def test_disabling_uniqueness_leaves_ids_alone(self, tmp_path):
        path = self._library(tmp_path, ["A", "A"])
        records = list(iter_library(path, unique_ids=False))
        assert [r.ligand_id for r in records] == ["A", "A"]


class TestSmilesReading:
    def test_reads_smiles_and_ids(self, tmp_path):
        path = _write(tmp_path, "l.smi", "CCO ethanol\nCCC propane\n")
        records = list(iter_library(path))

        assert [(r.ligand_id, r.block) for r in records] == [
            ("ethanol", "CCO"),
            ("propane", "CCC"),
        ]

    def test_skips_a_header_row(self, tmp_path):
        path = _write(tmp_path, "l.smi", "smiles name\nCCO ethanol\n")
        assert [r.ligand_id for r in iter_library(path)] == ["ethanol"]

    def test_unnamed_molecules_get_positional_ids(self, tmp_path):
        path = _write(tmp_path, "l.smi", "CCO\nCCC\n")
        assert [r.ligand_id for r in iter_library(path)] == ["mol1", "mol2"]

    def test_blank_and_comment_lines_are_ignored(self, tmp_path):
        path = _write(tmp_path, "l.smi", "# a comment\n\nCCO ethanol\n\n")
        assert [r.ligand_id for r in iter_library(path)] == ["ethanol"]


class TestMol2Reading:
    def test_splits_on_molecule_records(self, tmp_path):
        text = (
            "@<TRIPOS>MOLECULE\nfirst\n 1 0\n@<TRIPOS>ATOM\n"
            "@<TRIPOS>MOLECULE\nsecond\n 1 0\n@<TRIPOS>ATOM\n"
        )
        records = list(iter_library(_write(tmp_path, "l.mol2", text)))
        assert [r.ligand_id for r in records] == ["first", "second"]


class TestCounting:
    def test_counts_sdf_without_parsing(self, tmp_path):
        assert count_records(_write(tmp_path, "l.sdf", SDF_TEXT)) == 2

    def test_counts_an_unterminated_final_record(self, tmp_path):
        truncated = SDF_TEXT.rstrip()[: -len("$$$$")]
        assert count_records(_write(tmp_path, "l.sdf", truncated)) == 2

    def test_counts_smiles(self, tmp_path):
        path = _write(tmp_path, "l.smi", "CCO a\nCCC b\nCCCC c\n")
        assert count_records(path) == 3

    def test_count_agrees_with_iteration(self, tmp_path):
        path = _write(tmp_path, "l.sdf", SDF_TEXT * 3)
        assert count_records(path) == len(list(iter_library(path)))


class TestSurvey:
    def test_separates_records_from_distinct_compounds(self, tmp_path):
        blocks = []
        for cid in ["A", "A", "A", "B", "C", "C"]:
            blocks.append(
                "\n  test\n\n  1  0  0  0  0  0            999 V2000\n"
                "    0.0000    0.0000    0.0000 C   0  0\nM  END\n"
                f"> <COMPOUND_ID>\n{cid}\n\n$$$$\n"
            )
        path = _write(tmp_path, "l.sdf", "".join(blocks))

        info = survey(path)
        assert info["records"] == 6
        assert info["distinct_compounds"] == 3
        assert info["compounds_with_variants"] == 2
        assert info["max_variants"] == 3
