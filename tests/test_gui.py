"""Tests for the GUI's model and watcher.

The window itself is mostly layout and is left to manual inspection. What is
tested here is the part that has to hold up under a real screen: a table model
fed 47,000 rows, and a watcher that follows a journal being written underneath
it.

Everything runs offscreen (see conftest), so these execute on a headless runner.
"""

from __future__ import annotations

import pytest
from PySide6.QtCore import Qt

from drydock.core.rundir import LigandResult, PoseMode, RunDir
from drydock.core.synthetic import synthesize_run
from drydock.gui.model import ResultsTableModel
from drydock.gui.watcher import MAX_RECORDS_PER_POLL, RunWatcher


def _result(ligand_id: str, affinity: float | None = -7.5, elapsed: float = 1.0) -> LigandResult:
    if affinity is None:
        return LigandResult(ligand_id=ligand_id, status="failed", error="no poses", seed=1)
    return LigandResult(
        ligand_id=ligand_id,
        status="ok",
        seed=1,
        elapsed_s=elapsed,
        modes=(PoseMode(1, affinity), PoseMode(2, affinity + 0.4, 2.1, 3.3)),
    )


class TestResultsModel:
    def test_shape_matches_contents(self, qapp):
        model = ResultsTableModel()
        model.reset([_result("A"), _result("B")])

        assert model.rowCount() == 2
        assert model.columnCount() == len(ResultsTableModel.COLUMNS)

    def test_affinity_is_formatted_to_one_decimal(self, qapp):
        model = ResultsTableModel()
        model.reset([_result("A", -7.53)])

        index = model.index(0, 2)
        assert model.data(index, Qt.ItemDataRole.DisplayRole) == "-7.5"

    def test_failed_ligand_shows_a_dash_not_a_crash(self, qapp):
        model = ResultsTableModel()
        model.reset([_result("BAD", None)])

        assert model.data(model.index(0, 2), Qt.ItemDataRole.DisplayRole) == "-"
        assert model.data(model.index(0, 4), Qt.ItemDataRole.DisplayRole) == "failed"

    def test_error_surfaces_as_a_tooltip(self, qapp):
        model = ResultsTableModel()
        model.reset([_result("BAD", None)])

        assert "no poses" in model.data(model.index(0, 1), Qt.ItemDataRole.ToolTipRole)

    def test_sorts_by_affinity_numerically_not_lexically(self, qapp):
        """'-10.2' sorts before '-9.1' as text, and after it as a number."""
        model = ResultsTableModel()
        model.reset([_result("A", -9.1), _result("B", -10.2), _result("C", -3.0)])

        model.sort(2, Qt.SortOrder.AscendingOrder)
        assert [r.ligand_id for r in model.records] == ["B", "A", "C"]

    def test_failed_ligands_sort_to_the_end_either_way(self, qapp):
        """A missing affinity must never outrank a real one."""
        model = ResultsTableModel()
        model.reset([_result("OK1", -7.0), _result("BAD", None), _result("OK2", -9.0)])

        model.sort(2, Qt.SortOrder.AscendingOrder)
        assert model.records[-1].ligand_id == "BAD"

        model.sort(2, Qt.SortOrder.DescendingOrder)
        assert model.records[0].ligand_id == "BAD"
        assert [r.ligand_id for r in model.records[1:]] == ["OK1", "OK2"]

    def test_rank_column_renumbers_after_sorting(self, qapp):
        model = ResultsTableModel()
        model.reset([_result("A", -5.0), _result("B", -8.0)])
        model.sort(2, Qt.SortOrder.AscendingOrder)

        assert model.data(model.index(0, 0), Qt.ItemDataRole.DisplayRole) == "1"
        assert model.data(model.index(0, 1), Qt.ItemDataRole.DisplayRole) == "B"

    def test_append_keeps_the_active_sort(self, qapp):
        model = ResultsTableModel()
        model.reset([_result("A", -7.0)])
        model.sort(2, Qt.SortOrder.AscendingOrder)

        model.append([_result("BEST", -12.0), _result("WORST", -2.0)])

        assert [r.ligand_id for r in model.records] == ["BEST", "A", "WORST"]

    def test_append_of_nothing_is_a_no_op(self, qapp):
        model = ResultsTableModel()
        model.reset([_result("A")])
        model.append([])
        assert model.rowCount() == 1

    def test_reset_replaces_rather_than_accumulates(self, qapp):
        model = ResultsTableModel()
        model.reset([_result("A"), _result("B")])
        model.reset([_result("C")])

        assert model.rowCount() == 1
        assert model.records[0].ligand_id == "C"

    def test_out_of_range_access_returns_none(self, qapp):
        model = ResultsTableModel()
        model.reset([_result("A")])
        assert model.record_at(99) is None
        assert model.data(model.index(99, 0)) is None

    @pytest.mark.parametrize("n", [47_451])
    def test_holds_a_full_library_without_materialising_rows(self, qapp, n):
        """The whole point of QAbstractTableModel over QTableWidget.

        Qt asks for cells it is about to paint, so the cost of a large result set
        is the list of Python objects and nothing more.
        """
        model = ResultsTableModel()
        model.reset([_result(f"CMNPD{i}", -5.0 - (i % 50) / 10) for i in range(n)])

        assert model.rowCount() == n
        # Cells anywhere in the range must be reachable without touching the rest.
        assert model.data(model.index(0, 1), Qt.ItemDataRole.DisplayRole).startswith("CMNPD")
        assert model.data(model.index(n - 1, 1), Qt.ItemDataRole.DisplayRole).startswith("CMNPD")


class TestRunWatcher:
    def test_attach_to_missing_directory_reports_and_returns_false(self, qapp, tmp_path):
        watcher = RunWatcher()
        errors: list[str] = []
        watcher.error.connect(errors.append)

        assert watcher.attach(str(tmp_path / "nope")) is False
        assert errors and "No run directory" in errors[0]

    def test_attach_emits_status(self, qapp, tmp_path):
        synthesize_run(tmp_path / "run", n_ligands=20, seed=1)
        watcher = RunWatcher()
        seen = []
        watcher.statusChanged.connect(seen.append)

        assert watcher.attach(str(tmp_path / "run")) is True
        assert seen and seen[-1].total == 20
        watcher.detach()

    def test_load_all_returns_the_whole_journal(self, qapp, tmp_path):
        synthesize_run(tmp_path / "run", n_ligands=40, seed=2)
        watcher = RunWatcher()
        watcher.attach(str(tmp_path / "run"))

        assert len(watcher.load_all()) == 40
        watcher.detach()

    def test_load_all_leaves_no_records_to_replay(self, qapp, tmp_path):
        """Attaching must not deliver the same record twice."""
        synthesize_run(tmp_path / "run", n_ligands=30, seed=3)
        watcher = RunWatcher()
        watcher.attach(str(tmp_path / "run"))
        watcher.load_all()

        delivered: list[LigandResult] = []
        watcher.recordsAdded.connect(delivered.extend)
        watcher.poll()

        assert delivered == []
        watcher.detach()

    def test_poll_picks_up_records_written_after_attach(self, qapp, tmp_path):
        run = synthesize_run(tmp_path / "run", n_ligands=5, seed=4)
        watcher = RunWatcher()
        watcher.attach(str(run.path))
        watcher.load_all()

        delivered: list[LigandResult] = []
        watcher.recordsAdded.connect(delivered.extend)

        RunDir(run.path).append(_result("LATE", -9.9))
        watcher.poll()

        assert [r.ligand_id for r in delivered] == ["LATE"]
        watcher.detach()

    def test_large_batches_are_released_over_several_polls(self, qapp, tmp_path):
        """A run that advanced while the GUI was closed must not stall it."""
        run = RunDir(tmp_path / "run").create()
        n = MAX_RECORDS_PER_POLL + 500
        for i in range(n):
            run.append(_result(f"LIG{i}"))

        watcher = RunWatcher()
        batches: list[int] = []
        watcher.recordsAdded.connect(lambda recs: batches.append(len(recs)))

        # attach() polls once itself, which is the first capped batch.
        watcher.attach(str(run.path))
        assert batches == [MAX_RECORDS_PER_POLL], "no single poll may exceed the cap"

        watcher.poll()
        assert batches == [MAX_RECORDS_PER_POLL, 500]
        assert sum(batches) == n, "every record must arrive eventually"
        watcher.detach()

    def test_run_finished_fires_once_on_the_transition(self, qapp, tmp_path):
        synthesize_run(tmp_path / "run", n_ligands=10, seed=5)
        watcher = RunWatcher()
        finished = []
        watcher.runFinished.connect(finished.append)

        watcher.attach(str(tmp_path / "run"))
        watcher.poll()
        watcher.poll()

        assert len(finished) == 1
        assert finished[0].state == "finished"
        watcher.detach()

    def test_incomplete_run_does_not_report_finished(self, qapp, tmp_path):
        synthesize_run(tmp_path / "run", n_ligands=50, seed=6, complete=False)
        watcher = RunWatcher()
        finished = []
        watcher.runFinished.connect(finished.append)

        watcher.attach(str(tmp_path / "run"))
        assert finished == []
        watcher.detach()

    def test_detach_stops_delivering(self, qapp, tmp_path):
        run = synthesize_run(tmp_path / "run", n_ligands=5, seed=7)
        watcher = RunWatcher()
        watcher.attach(str(run.path))
        watcher.detach()

        delivered: list[LigandResult] = []
        watcher.recordsAdded.connect(delivered.extend)
        RunDir(run.path).append(_result("AFTER"))
        watcher.poll()

        assert delivered == []
        assert watcher.run is None

    def test_survives_a_missing_status_file(self, qapp, tmp_path):
        """status.json is a cache; losing it must not stop the watcher."""
        run = synthesize_run(tmp_path / "run", n_ligands=10, seed=8)
        run.status_file.unlink()

        watcher = RunWatcher()
        seen = []
        watcher.statusChanged.connect(seen.append)
        watcher.attach(str(run.path))

        assert seen, "status must be rebuilt from the journal"
        assert seen[-1].completed + seen[-1].failed == 10
        watcher.detach()


class TestMainWindow:
    """Window-level tests.

    These exist because a model that sorts correctly in isolation can still be
    re-sorted by the view during construction. ``setSortingEnabled(True)`` applies
    the header's default indicator, column 0, *descending*, which presented
    the worst binders at the top of a screening result. The model unit tests all
    passed throughout.
    """

    def test_opens_with_best_affinity_first(self, qapp, tmp_path):
        from drydock.gui.app import MainWindow

        run = synthesize_run(tmp_path / "run", n_ligands=300, seed=11, failure_rate=0.05)
        window = MainWindow()
        window.attach_run(str(run.path))

        rows = window._model.records
        scored = [r.best_affinity for r in rows if r.best_affinity is not None]
        assert scored, "synthetic run should contain successful ligands"
        assert scored == sorted(scored), "best (most negative) affinity must lead"
        assert rows[0].best_affinity == min(scored)

    def test_failed_ligands_do_not_occupy_the_top(self, qapp, tmp_path):
        from drydock.gui.app import MainWindow

        run = synthesize_run(tmp_path / "run", n_ligands=300, seed=12, failure_rate=0.3)
        window = MainWindow()
        window.attach_run(str(run.path))

        assert window._model.records[0].status == "ok"

    def test_header_indicator_matches_the_model_sort(self, qapp, tmp_path):
        """A header arrow pointing the other way to the data is its own bug."""
        from drydock.gui.app import AFFINITY_COLUMN, MainWindow

        run = synthesize_run(tmp_path / "run", n_ligands=50, seed=13)
        window = MainWindow()
        window.attach_run(str(run.path))

        header = window._table.horizontalHeader()
        assert header.sortIndicatorSection() == AFFINITY_COLUMN
        assert header.sortIndicatorOrder() == Qt.SortOrder.AscendingOrder

    def test_attaching_to_a_second_run_replaces_the_first(self, qapp, tmp_path):
        from drydock.gui.app import MainWindow

        first = synthesize_run(tmp_path / "a", n_ligands=40, seed=14)
        second = synthesize_run(tmp_path / "b", n_ligands=10, seed=15)

        window = MainWindow()
        window.attach_run(str(first.path))
        window.attach_run(str(second.path))

        assert window._model.rowCount() == 10

    def test_detach_leaves_the_window_usable(self, qapp, tmp_path):
        from drydock.gui.app import MainWindow

        run = synthesize_run(tmp_path / "run", n_ligands=20, seed=16)
        window = MainWindow()
        window.attach_run(str(run.path))
        window.detach_run()

        assert window._watcher.run is None


class TestWatcherModelIntegration:
    def test_records_flow_from_journal_into_the_model(self, qapp, tmp_path):
        run = synthesize_run(tmp_path / "run", n_ligands=100, seed=9, failure_rate=0.1)

        model = ResultsTableModel()
        watcher = RunWatcher()
        watcher.recordsAdded.connect(model.append)
        watcher.attach(str(run.path))
        watcher.poll()

        assert model.rowCount() == 100
        affinities = [
            r.best_affinity for r in model.records if r.best_affinity is not None
        ]
        assert affinities == sorted(affinities), "model should present best-first"
        watcher.detach()
