"""Tests for results assembly.

The behaviour that matters most here is stereoisomer grouping. A library that
enumerates every stereoisomer under one identifier -- CMNPD gives CMNPD22318 64
of them -- will otherwise fill a top-100 hit list with variants of a handful of
compounds, which looks like a strong result and is not one.
"""

from __future__ import annotations

import csv

import pytest

from drydock.core.results import ALL_MODES_COLUMNS, RESULT_COLUMNS, collate, write_results
from drydock.core.rundir import LigandResult, PoseMode, RunDir


def _result(ligand_id: str, affinity: float | None, seed: int = 42) -> LigandResult:
    if affinity is None:
        return LigandResult(ligand_id=ligand_id, status="failed", seed=seed, error="no poses")
    return LigandResult(
        ligand_id=ligand_id,
        status="ok",
        seed=seed,
        elapsed_s=1.0,
        modes=(
            PoseMode(1, affinity, 0.0, 0.0),
            PoseMode(2, affinity + 0.5, 2.1, 3.4),
        ),
    )


def _manifest(**entries: tuple[str, int]) -> dict[str, dict[str, str]]:
    """Build a manifest keyed by ligand_id from (compound_id, heavy_atoms)."""
    return {
        ligand_id: {
            "ligand_id": ligand_id,
            "compound_id": compound_id,
            "heavy_atoms": str(heavy),
            "mw": "350.0",
            "clogp": "3.2",
            "formula": "C20H30O2",
            "smiles": "CCO",
        }
        for ligand_id, (compound_id, heavy) in entries.items()
    }


class TestCollate:
    def test_ranks_best_first(self):
        rows = collate(
            [_result("A", -7.0), _result("B", -9.5), _result("C", -5.0)],
            _manifest(A=("A", 20), B=("B", 20), C=("C", 20)),
        )
        assert [r.ligand_id for r in rows] == ["B", "A", "C"]

    def test_failed_ligands_sort_last(self):
        rows = collate(
            [_result("BAD", None), _result("GOOD", -8.0)],
            _manifest(BAD=("BAD", 20), GOOD=("GOOD", 20)),
        )
        assert rows[0].ligand_id == "GOOD"
        assert rows[-1].best_affinity is None

    def test_computes_ligand_efficiency(self):
        rows = collate([_result("A", -10.0)], _manifest(A=("A", 25)))
        assert rows[0].ligand_efficiency == pytest.approx(-0.4)

    def test_ligand_efficiency_is_none_without_a_manifest(self):
        """Heavy-atom count comes from preparation; without it, LE is unknowable."""
        rows = collate([_result("A", -10.0)], {})
        assert rows[0].ligand_efficiency is None

    def test_missing_manifest_entry_still_produces_a_row(self):
        rows = collate([_result("ORPHAN", -8.0)], {})
        assert rows[0].compound_id == "ORPHAN"
        assert rows[0].best_affinity == -8.0


class TestStereoisomerGrouping:
    def test_keeps_only_the_best_variant(self):
        records = [
            _result("CPD", -7.0),
            _result("CPD_2", -9.5),
            _result("CPD_3", -6.0),
        ]
        manifest = _manifest(CPD=("CPD", 20), CPD_2=("CPD", 20), CPD_3=("CPD", 20))

        rows = collate(records, manifest, group_stereoisomers=True)
        assert len(rows) == 1
        assert rows[0].best_affinity == -9.5
        assert rows[0].ligand_id == "CPD_2", "the winning variant must be identifiable"
        assert rows[0].n_variants == 3

    def test_flat_mode_keeps_every_record(self):
        records = [_result("CPD", -7.0), _result("CPD_2", -9.5)]
        manifest = _manifest(CPD=("CPD", 20), CPD_2=("CPD", 20))

        rows = collate(records, manifest, group_stereoisomers=False)
        assert len(rows) == 2

    def test_grouping_does_not_let_one_compound_dominate(self):
        """The behaviour this exists for."""
        records = [_result(f"BIG_{i}", -9.0 - i * 0.01) for i in range(64)]
        records.append(_result("OTHER", -8.9))
        manifest = _manifest(
            **{f"BIG_{i}": ("BIG", 30) for i in range(64)}, OTHER=("OTHER", 20)
        )

        grouped = collate(records, manifest, group_stereoisomers=True)
        assert len(grouped) == 2, "64 variants collapse to one compound"
        assert grouped[0].compound_id == "BIG"
        assert grouped[0].n_variants == 64

        flat = collate(records, manifest, group_stereoisomers=False)
        assert len(flat) == 65
        assert all(r.compound_id == "BIG" for r in flat[:10]), (
            "ungrouped, one compound occupies the whole head of the list"
        )

    def test_a_failed_variant_does_not_beat_a_scored_one(self):
        records = [_result("CPD", None), _result("CPD_2", -7.0)]
        manifest = _manifest(CPD=("CPD", 20), CPD_2=("CPD", 20))

        rows = collate(records, manifest, group_stereoisomers=True)
        assert rows[0].best_affinity == -7.0
        assert rows[0].ligand_id == "CPD_2"


class TestWriteResults:
    def test_writes_both_files(self, tmp_path):
        run = RunDir(tmp_path / "run").create()
        run.append(_result("A", -8.0))
        run.append(_result("B", -6.0))

        results, all_modes, n = write_results(run.path)

        assert results.exists() and all_modes.exists()
        assert n == 2

    def test_results_columns_are_stable(self, tmp_path):
        run = RunDir(tmp_path / "run").create()
        run.append(_result("A", -8.0))
        results, _, _ = write_results(run.path)

        with open(results, newline="", encoding="utf-8") as fh:
            assert tuple(csv.DictReader(fh).fieldnames) == RESULT_COLUMNS

    def test_all_modes_matches_the_padel_schema(self, tmp_path):
        """Column names and order are what old downstream scripts expect."""
        run = RunDir(tmp_path / "run").create()
        run.append(_result("A", -8.0))
        _, all_modes, _ = write_results(run.path)

        with open(all_modes, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            assert tuple(reader.fieldnames) == ALL_MODES_COLUMNS
            rows = list(reader)

        assert len(rows) == 2, "one row per mode"
        assert rows[0]["Ligand"] == "A"
        assert rows[0]["Mode"] == "1"
        assert rows[0]["Affinity"] == "-8.0"
        assert rows[0]["RandomSeed"] == "42"

    def test_ranks_are_sequential_from_one(self, tmp_path):
        run = RunDir(tmp_path / "run").create()
        for i, affinity in enumerate([-5.0, -9.0, -7.0]):
            run.append(_result(f"L{i}", affinity))

        results, _, _ = write_results(run.path)
        with open(results, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))

        assert [r["rank"] for r in rows] == ["1", "2", "3"]
        assert rows[0]["best_affinity"] == "-9.0"

    def test_joins_descriptors_from_a_manifest(self, tmp_path):
        run = RunDir(tmp_path / "run").create()
        run.append(_result("LIG", -10.0))

        manifest_path = tmp_path / "manifest.csv"
        with open(manifest_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=["ligand_id", "compound_id", "heavy_atoms", "mw", "formula"]
            )
            writer.writeheader()
            writer.writerow(
                {
                    "ligand_id": "LIG",
                    "compound_id": "CPD",
                    "heavy_atoms": "25",
                    "mw": "340.5",
                    "formula": "C20H28O4",
                }
            )

        results, _, _ = write_results(run.path, manifest_path)
        with open(results, newline="", encoding="utf-8") as fh:
            row = next(csv.DictReader(fh))

        assert row["compound_id"] == "CPD"
        assert row["mw"] == "340.5"
        assert row["ligand_efficiency"] == "-0.4"

    def test_can_report_on_a_run_in_progress(self, tmp_path):
        """Reading the journal directly means no need to wait for the run."""
        run = RunDir(tmp_path / "run").create()
        run.append(_result("A", -8.0))

        _, _, n = write_results(run.path)
        assert n == 1

        run.append(_result("B", -9.0))
        _, _, n = write_results(run.path)
        assert n == 2

    def test_empty_run_writes_headers_only(self, tmp_path):
        run = RunDir(tmp_path / "run").create()
        results, all_modes, n = write_results(run.path)

        assert n == 0
        with open(results, newline="", encoding="utf-8") as fh:
            assert list(csv.DictReader(fh)) == []
