"""Tests for provenance recording and the self-test.

The self-test's own comparison logic is tested here rather than by running the
pipeline, which takes ~18 seconds and is exercised by ``drydock selftest`` and by
CI. What matters is that the comparison *fails* when it should: a check that can
only pass verifies nothing, and a reference comparison that silently tolerates
drift is worse than no comparison at all.
"""

from __future__ import annotations

import json

import pytest

from drydock.core import provenance
from drydock.core.selftest import (
    AFFINITY_TOLERANCE,
    LIGANDS,
    RECEPTOR,
    REFERENCE,
    SelftestResult,
    _compare_to_reference,
    load_reference,
)


class TestBundledData:
    def test_receptor_is_present_and_small(self):
        """Shipped in the package, so it must stay small enough to be reasonable."""
        assert RECEPTOR.exists()
        assert RECEPTOR.stat().st_size < 200_000

    def test_receptor_uses_canonical_metal_casing(self):
        """'ZN' would make the bundled case fail on every installation."""
        from drydock.core.receptor import inspect

        report = inspect(RECEPTOR)
        assert "Zn" in report.atom_types
        assert "ZN" not in report.atom_types

    def test_ligands_are_present(self):
        assert LIGANDS.exists()
        lines = [ln for ln in LIGANDS.read_text().splitlines() if ln.strip()]
        assert len(lines) >= 3

    def test_reference_is_recorded(self):
        reference = load_reference()
        assert reference is not None, "run 'drydock dev record-reference'"
        assert reference["affinities"]

    def test_reference_records_the_settings_it_was_measured_with(self):
        """A reference only means something alongside how it was produced."""
        reference = load_reference()
        settings = reference["settings"]
        assert settings["seed"] is not None
        assert settings["exhaustiveness"] >= 1
        assert len(settings["box_center"]) == 3

    def test_reference_records_engine_versions(self):
        reference = load_reference()
        assert reference["packages"]["vina"]
        assert reference["packages"]["rdkit"]


class TestReferenceComparison:
    """The comparison must fail when results drift. Otherwise it is decoration."""

    def _reference(self, **affinities: float) -> dict:
        return {"affinities": affinities}

    def test_passes_on_exact_agreement(self):
        result = SelftestResult()
        result.affinities = {"a": -7.5, "b": -5.0}
        _compare_to_reference(result, self._reference(a=-7.5, b=-5.0))

        assert result.passed

    def test_passes_within_tolerance(self):
        result = SelftestResult()
        result.affinities = {"a": -7.5 + AFFINITY_TOLERANCE * 0.9}
        _compare_to_reference(result, self._reference(a=-7.5))

        assert result.passed

    def test_fails_outside_tolerance(self):
        result = SelftestResult()
        result.affinities = {"a": -7.5 + AFFINITY_TOLERANCE * 2}
        _compare_to_reference(result, self._reference(a=-7.5))

        assert not result.passed
        assert any("within" in c.name for c in result.failures)

    def test_failure_names_the_offending_ligand_and_both_values(self):
        result = SelftestResult()
        result.affinities = {"good": -5.0, "drifted": -1.0}
        _compare_to_reference(result, self._reference(good=-5.0, drifted=-9.0))

        detail = " ".join(c.detail for c in result.failures)
        assert "drifted" in detail
        assert "-9.0" in detail and "-1.0" in detail
        assert "good" not in detail

    def test_fails_when_a_ligand_is_missing_entirely(self):
        """A ligand that stopped docking must not pass by being absent."""
        result = SelftestResult()
        result.affinities = {"a": -7.5}
        _compare_to_reference(result, self._reference(a=-7.5, b=-5.0))

        assert not result.passed
        assert any("b" in c.detail for c in result.failures)

    def test_fails_on_an_empty_reference(self):
        result = SelftestResult()
        result.affinities = {"a": -7.5}
        _compare_to_reference(result, {"affinities": {}})

        assert not result.passed

    @pytest.mark.parametrize("delta", [-2.0, -1.0, 1.0, 2.0])
    def test_drift_fails_in_either_direction(self, delta):
        """A score improving unexpectedly is as much a change as one worsening."""
        result = SelftestResult()
        result.affinities = {"a": -7.5 + delta}
        _compare_to_reference(result, self._reference(a=-7.5))

        assert not result.passed


class TestSelftestResult:
    def test_passed_requires_every_check(self):
        result = SelftestResult()
        result.add("one", True)
        result.add("two", True)
        assert result.passed

        result.add("three", False)
        assert not result.passed

    def test_failures_lists_only_failures(self):
        result = SelftestResult()
        result.add("ok", True)
        result.add("bad", False, "why")
        assert [c.name for c in result.failures] == ["bad"]


class TestProvenance:
    def test_hashes_a_file(self, tmp_path):
        path = tmp_path / "f.txt"
        path.write_text("hello", encoding="utf-8")
        digest = provenance.file_sha256(path)

        assert digest == (
            "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824"
        )

    def test_missing_file_hashes_to_none(self, tmp_path):
        assert provenance.file_sha256(tmp_path / "absent") is None

    def test_hash_detects_a_changed_receptor(self, tmp_path):
        """Why inputs are hashed rather than just named."""
        path = tmp_path / "rec.pdbqt"
        path.write_text("ATOM      1  N\n", encoding="utf-8")
        before = provenance.file_sha256(path)

        path.write_text("ATOM      1  C\n", encoding="utf-8")
        assert provenance.file_sha256(path) != before

    def test_directory_digest_summarises_without_hashing_every_file(self, tmp_path):
        for i in range(5):
            (tmp_path / f"lig{i}.pdbqt").write_text("x" * (i + 1), encoding="utf-8")

        digest = provenance.directory_digest(tmp_path)
        assert digest["n_files"] == 5
        assert digest["total_bytes"] == 1 + 2 + 3 + 4 + 5
        assert len(digest["listing_sha256"]) == 64

    def test_directory_digest_changes_when_a_file_is_added(self, tmp_path):
        (tmp_path / "a.pdbqt").write_text("a", encoding="utf-8")
        before = provenance.directory_digest(tmp_path)["listing_sha256"]

        (tmp_path / "b.pdbqt").write_text("b", encoding="utf-8")
        assert provenance.directory_digest(tmp_path)["listing_sha256"] != before

    def test_directory_digest_changes_when_a_file_is_resized(self, tmp_path):
        path = tmp_path / "a.pdbqt"
        path.write_text("a", encoding="utf-8")
        before = provenance.directory_digest(tmp_path)["listing_sha256"]

        path.write_text("aa", encoding="utf-8")
        assert provenance.directory_digest(tmp_path)["listing_sha256"] != before

    def test_missing_directory_is_reported_not_raised(self, tmp_path):
        assert provenance.directory_digest(tmp_path / "absent")["exists"] is False

    def test_records_tracked_package_versions(self):
        versions = provenance.package_versions()
        assert set(provenance.TRACKED_PACKAGES) <= set(versions)
        assert versions["vina"] != "not installed"
        assert versions["rdkit"] != "not installed"

    def test_records_external_engine_versions(self):
        """autogrid4's version is invisible to Python packaging but sets results."""
        versions = provenance.engine_versions()
        assert "4.2.6" in versions["autogrid4"]

    def test_build_produces_a_serialisable_record(self, tmp_path):
        from drydock.core.box import Box
        from drydock.engines.base import DockConfig

        receptor = tmp_path / "rec.pdbqt"
        receptor.write_text("ATOM      1  N\n", encoding="utf-8")
        ligands = tmp_path / "ligs"
        ligands.mkdir()

        config = DockConfig(
            receptor=str(receptor), box=Box((0.0, 0.0, 0.0), (20.0, 20.0, 20.0)), seed=7
        )
        record = provenance.build(receptor=receptor, ligand_dir=ligands, config=config)

        json.dumps(record)  # must not raise
        assert record["receptor"]["sha256"]
        assert record["reproducibility"]["global_seed"] == 7
        assert record["reproducibility"]["threads_per_job"] == 1
        assert record["config"]["box"]["size"] == [20.0, 20.0, 20.0]

    def test_build_explains_the_single_thread_choice(self):
        """The claim a reader most needs to evaluate should not be implicit."""
        from drydock.core.box import Box
        from drydock.engines.base import DockConfig

        config = DockConfig(receptor="r.pdbqt", box=Box((0.0, 0.0, 0.0), (20.0, 20.0, 20.0)))
        record = provenance.build(receptor="r.pdbqt", ligand_dir=".", config=config)

        assert "thread" in record["reproducibility"]["note"].lower()


def test_reference_file_is_valid_json():
    assert json.loads(REFERENCE.read_text())["affinities"]
