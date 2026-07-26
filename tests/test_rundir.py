"""Tests for the run-directory contract.

Weighted towards the failure modes the design exists to survive -- truncated
journals, stale caches, resume after a kill -- rather than the happy path, which
is exercised incidentally by everything else.
"""

from __future__ import annotations

import json

import pytest

from drydock.core.rundir import LigandResult, PoseMode, RunDir, RunStatus
from drydock.core.synthetic import synthesize_run


@pytest.fixture
def run(tmp_path):
    return RunDir(tmp_path / "run").create()


def _result(ligand_id: str, affinity: float = -7.5, status: str = "ok") -> LigandResult:
    modes = (
        (
            PoseMode(1, affinity, 0.0, 0.0),
            PoseMode(2, affinity + 0.3, 1.9, 3.2),
        )
        if status == "ok"
        else ()
    )
    return LigandResult(ligand_id=ligand_id, status=status, seed=42, elapsed_s=1.0, modes=modes)


class TestJournal:
    def test_roundtrip_preserves_records(self, run):
        run.append(_result("LIG1", -7.5))
        run.append(_result("LIG2", -6.1))

        got = list(run.read_journal())
        assert [r.ligand_id for r in got] == ["LIG1", "LIG2"]
        assert got[0].best_affinity == -7.5
        assert got[0].modes[1].rmsd_ub == 3.2

    def test_best_affinity_is_the_minimum_not_the_first(self, run):
        """Engines normally emit sorted modes, but the contract must not rely on it."""
        result = LigandResult(
            ligand_id="LIG",
            status="ok",
            modes=(PoseMode(1, -5.0), PoseMode(2, -8.2), PoseMode(3, -6.0)),
        )
        assert result.best_affinity == -8.2

    def test_no_modes_means_no_affinity(self):
        assert LigandResult(ligand_id="LIG", status="failed").best_affinity is None

    def test_truncated_final_line_is_skipped(self, run):
        """The signature of a killed writer: everything before it must survive."""
        run.append(_result("LIG1"))
        run.append(_result("LIG2"))
        with open(run.journal_file, "a", encoding="utf-8") as fh:
            fh.write('{"ligand_id": "LIG3", "status": "o')

        got = list(run.read_journal())
        assert [r.ligand_id for r in got] == ["LIG1", "LIG2"]

    def test_garbage_lines_do_not_abort_the_read(self, run):
        run.append(_result("LIG1"))
        with open(run.journal_file, "a", encoding="utf-8") as fh:
            fh.write("not json at all\n")
            fh.write('{"missing": "required fields"}\n')
        run.append(_result("LIG2"))

        assert [r.ligand_id for r in run.read_journal()] == ["LIG1", "LIG2"]

    def test_completed_ids_includes_failures(self, run):
        """A ligand that reliably kills the engine must not be retried forever."""
        run.append(_result("OK", status="ok"))
        run.append(_result("BAD", status="failed"))

        assert run.completed_ids() == {"OK", "BAD"}

    def test_missing_journal_reads_as_empty(self, tmp_path):
        assert list(RunDir(tmp_path / "absent").read_journal()) == []


class TestTailJournal:
    def test_returns_only_new_records(self, run):
        run.append(_result("LIG1"))
        first, offset = run.tail_journal(0)
        assert [r.ligand_id for r in first] == ["LIG1"]

        run.append(_result("LIG2"))
        second, offset2 = run.tail_journal(offset)
        assert [r.ligand_id for r in second] == ["LIG2"]
        assert offset2 > offset

        assert run.tail_journal(offset2) == ([], offset2)

    def test_partial_line_is_left_for_the_next_poll(self, run):
        """A watcher must never see a record the writer has not finished."""
        run.append(_result("LIG1"))
        _, offset = run.tail_journal(0)

        with open(run.journal_file, "a", encoding="utf-8") as fh:
            fh.write('{"ligand_id": "LIG2", "sta')
        records, new_offset = run.tail_journal(offset)
        assert records == []
        assert new_offset == offset, "offset must not advance past an incomplete record"

        with open(run.journal_file, "a", encoding="utf-8") as fh:
            fh.write('tus": "ok", "modes": [{"mode": 1, "affinity": -7.0}]}\n')
        records, _ = run.tail_journal(new_offset)
        assert [r.ligand_id for r in records] == ["LIG2"]

    def test_restarts_when_journal_shrinks(self, run):
        """A replaced journal must not be read from a now-meaningless offset."""
        run.append(_result("LIG1"))
        run.append(_result("LIG2"))
        _, offset = run.tail_journal(0)

        run.journal_file.unlink()
        run.append(_result("FRESH"))

        records, _ = run.tail_journal(offset)
        assert [r.ligand_id for r in records] == ["FRESH"]


class TestStatus:
    def test_roundtrip(self, run):
        run.write_status(RunStatus(state="running", total=100, completed=10, engine="vina"))
        got = run.read_status()
        assert (got.state, got.total, got.completed, got.engine) == ("running", 100, 10, "vina")
        assert got.updated_at is not None

    def test_corrupt_status_reads_as_none_not_error(self, run):
        """status.json is a cache; a bad one costs a rebuild, not a run."""
        run.status_file.write_text("{ this is not json", encoding="utf-8")
        assert run.read_status() is None

    def test_rebuild_recounts_from_the_journal(self, run):
        run.write_status(RunStatus(state="running", total=3, completed=999, engine="vina"))
        run.append(_result("A", status="ok"))
        run.append(_result("B", status="ok"))
        run.append(_result("C", status="failed"))

        rebuilt = run.rebuild_status()
        assert (rebuilt.completed, rebuilt.failed) == (2, 1)
        assert rebuilt.total == 3, "total and engine carry over from the cache"
        assert rebuilt.engine == "vina"

    def test_progress_arithmetic(self):
        status = RunStatus(state="running", total=100, completed=20, failed=5, skipped=5)
        assert status.done == 30
        assert status.remaining == 70
        assert status.fraction == pytest.approx(0.3)

    def test_rate_and_eta_are_none_before_progress(self):
        assert RunStatus(state="running", total=10).rate_per_s is None
        assert RunStatus(state="running", total=10).eta_s is None

    def test_eta_projects_from_observed_rate(self):
        import time as _time

        now = _time.time()
        status = RunStatus(
            state="running", total=100, completed=25, started_at=now - 50, updated_at=now
        )
        assert status.rate_per_s == pytest.approx(0.5, rel=1e-3)
        assert status.eta_s == pytest.approx(150.0, rel=1e-3)

    def test_atomic_write_leaves_no_temp_files(self, run):
        for _ in range(5):
            run.write_status(RunStatus(state="running", total=1))
        assert list(run.path.glob(".status.json.*")) == []


class TestProvenance:
    def test_roundtrip_stamps_schema_version(self, run):
        run.write_provenance({"engine": "vina", "global_seed": 7})
        got = run.read_provenance()
        assert got["engine"] == "vina"
        assert got["schema_version"] == 1

    def test_missing_provenance_reads_as_none(self, run):
        assert run.read_provenance() is None


class TestSyntheticRuns:
    def test_generates_a_complete_run(self, tmp_path):
        run = synthesize_run(tmp_path / "run", n_ligands=50, seed=1, failure_rate=0.0)

        records = list(run.read_journal())
        assert len(records) == 50
        assert all(r.status == "ok" for r in records)
        assert run.read_status().state == "finished"
        assert all(len(r.modes) == 9 for r in records)

    def test_is_reproducible_for_a_given_seed(self, tmp_path):
        a = synthesize_run(tmp_path / "a", n_ligands=25, seed=7)
        b = synthesize_run(tmp_path / "b", n_ligands=25, seed=7)

        affinities_a = [r.best_affinity for r in a.read_journal()]
        affinities_b = [r.best_affinity for r in b.read_journal()]
        assert affinities_a == affinities_b

    def test_incomplete_run_looks_like_a_killed_run(self, tmp_path):
        run = synthesize_run(tmp_path / "run", n_ligands=100, seed=2, complete=False)

        status = run.read_status()
        assert status.state == "running"
        assert status.total == 100
        assert 0 < len(list(run.read_journal())) < 100

    def test_failures_are_recorded_without_modes(self, tmp_path):
        run = synthesize_run(tmp_path / "run", n_ligands=200, seed=3, failure_rate=0.5)

        failed = [r for r in run.read_journal() if r.status == "failed"]
        assert failed, "failure_rate=0.5 should produce failures"
        assert all(r.modes == () and r.error for r in failed)

    def test_modes_worsen_monotonically(self, tmp_path):
        """Mode 1 is the best pose; the contract and the CSV both assume it."""
        run = synthesize_run(tmp_path / "run", n_ligands=20, seed=4, failure_rate=0.0)

        for record in run.read_journal():
            affinities = [m.affinity for m in record.modes]
            assert affinities == sorted(affinities)
            assert record.best_affinity == affinities[0]


class TestResumeFlow:
    def test_resume_skips_everything_already_recorded(self, tmp_path):
        run = synthesize_run(tmp_path / "run", n_ligands=100, seed=5, complete=False)
        done = run.completed_ids()

        all_ids = {f"CMNPD{i + 1}" for i in range(100)}
        todo = all_ids - done

        assert todo, "an incomplete run must leave work to do"
        assert todo.isdisjoint(done)
        assert len(todo) + len(done) == 100

    def test_appending_after_resume_preserves_earlier_records(self, tmp_path):
        run = synthesize_run(tmp_path / "run", n_ligands=20, seed=6, complete=False)
        before = len(list(run.read_journal()))

        RunDir(run.path).append(_result("RESUMED"))

        records = list(run.read_journal())
        assert len(records) == before + 1
        assert records[-1].ligand_id == "RESUMED"

    def test_journal_is_valid_jsonl(self, tmp_path):
        """Downstream tooling should be able to treat it as plain JSON Lines."""
        run = synthesize_run(tmp_path / "run", n_ligands=30, seed=8)

        for line in run.journal_file.read_text(encoding="utf-8").splitlines():
            assert json.loads(line)["ligand_id"]
