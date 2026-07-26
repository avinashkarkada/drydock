"""Turning a journal into a hit list.

Screening produces identifiers and affinities. On its own that is not a result:
``CMNPD31204, -11.8`` tells you nothing about whether the compound is worth
chasing. This module joins those numbers back to the chemistry recorded during
preparation and ranks them.

Two files come out, because two different questions get asked of a screen.

``results.csv``
    One row per compound, ranked. What you read to decide what to test.

``results_all_modes.csv``
    Every pose of every ligand, in the exact schema PaDEL-ADV emitted
    (``Ligand,Mode,Affinity,RMSD_LB,RMSD_UB,RandomSeed``), so existing downstream
    analyses keep working unchanged.

Grouping stereoisomers
----------------------

Libraries that enumerate stereoisomers under one identifier -- CMNPD gives
``CMNPD22318`` 64 of them -- will otherwise fill a top-100 list with variants of
a handful of compounds. Results are therefore grouped by ``compound_id`` and
ranked on the best-scoring variant, with the number of variants tried recorded
alongside. ``--flat`` reports every record separately instead.
"""

from __future__ import annotations

import csv
import os
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from drydock.core.descriptors import ligand_efficiency
from drydock.core.rundir import LigandResult, RunDir

# Columns of results.csv, in order. Identifiers, then the result, then the
# chemistry needed to judge it.
RESULT_COLUMNS: tuple[str, ...] = (
    "rank",
    "compound_id",
    "ligand_id",
    "best_affinity",
    "ligand_efficiency",
    "mw",
    "clogp",
    "tpsa",
    "heavy_atoms",
    "rot_bonds",
    "hbd",
    "hba",
    "formal_charge",
    "torsions",
    "n_variants",
    "n_modes",
    "formula",
    "smiles",
    "status",
    "seed",
)

# PaDEL-ADV's exact schema. Column names and order are load-bearing: they are
# what downstream scripts written against the old tool expect to find.
ALL_MODES_COLUMNS: tuple[str, ...] = (
    "Ligand",
    "Mode",
    "Affinity",
    "RMSD_LB",
    "RMSD_UB",
    "RandomSeed",
)


@dataclass(slots=True)
class ResultRow:
    """One ranked compound."""

    compound_id: str
    ligand_id: str
    best_affinity: float | None
    status: str
    seed: int | None
    n_modes: int
    n_variants: int
    descriptors: dict[str, Any]

    @property
    def heavy_atoms(self) -> int | None:
        value = self.descriptors.get("heavy_atoms")
        try:
            return int(value) if value not in (None, "") else None
        except (TypeError, ValueError):
            return None

    @property
    def ligand_efficiency(self) -> float | None:
        return ligand_efficiency(self.best_affinity, self.heavy_atoms)

    def to_row(self, rank: int) -> dict[str, Any]:
        desc = self.descriptors
        return {
            "rank": rank,
            "compound_id": self.compound_id,
            "ligand_id": self.ligand_id,
            "best_affinity": self.best_affinity,
            "ligand_efficiency": self.ligand_efficiency,
            "mw": desc.get("mw"),
            "clogp": desc.get("clogp"),
            "tpsa": desc.get("tpsa"),
            "heavy_atoms": desc.get("heavy_atoms"),
            "rot_bonds": desc.get("rot_bonds"),
            "hbd": desc.get("hbd"),
            "hba": desc.get("hba"),
            "formal_charge": desc.get("formal_charge"),
            "torsions": desc.get("torsions"),
            "n_variants": self.n_variants,
            "n_modes": self.n_modes,
            "formula": desc.get("formula"),
            "smiles": desc.get("smiles"),
            "status": self.status,
            "seed": self.seed,
        }


def _sort_key(row: ResultRow) -> tuple[int, float]:
    """Rank by affinity, with unscored compounds last regardless of direction."""
    if row.best_affinity is None:
        return (1, 0.0)
    return (0, row.best_affinity)


def collate(
    records: Iterable[LigandResult],
    manifest: dict[str, dict[str, str]],
    group_stereoisomers: bool = True,
) -> list[ResultRow]:
    """Join docking records to the ligand manifest and rank them.

    Args:
        records: Journal records.
        manifest: Ligand manifest keyed by ``ligand_id``.
        group_stereoisomers: Collapse records sharing a ``compound_id``, keeping
            the best-scoring variant.

    Returns:
        Rows sorted best-first.
    """
    rows: list[ResultRow] = []
    for record in records:
        entry = manifest.get(record.ligand_id, {})
        rows.append(
            ResultRow(
                compound_id=entry.get("compound_id") or record.ligand_id,
                ligand_id=record.ligand_id,
                best_affinity=record.best_affinity,
                status=record.status,
                seed=record.seed,
                n_modes=len(record.modes),
                n_variants=1,
                descriptors=dict(entry),
            )
        )

    if group_stereoisomers:
        best: dict[str, ResultRow] = {}
        counts: dict[str, int] = {}
        for row in rows:
            counts[row.compound_id] = counts.get(row.compound_id, 0) + 1
            incumbent = best.get(row.compound_id)
            if incumbent is None or _sort_key(row) < _sort_key(incumbent):
                best[row.compound_id] = row
        rows = list(best.values())
        for row in rows:
            row.n_variants = counts[row.compound_id]

    rows.sort(key=_sort_key)
    return rows


def write_results(
    run_dir: str | os.PathLike[str],
    manifest_path: str | os.PathLike[str] | None = None,
    group_stereoisomers: bool = True,
) -> tuple[Path, Path, int]:
    """Write ``results.csv`` and ``results_all_modes.csv``.

    Reads the journal directly, so it is safe to call against a run in progress
    and produces a hit list of whatever has finished so far.

    Returns:
        The two paths written and the number of ranked rows.
    """
    run = RunDir(run_dir)
    records = list(run.read_journal())

    manifest: dict[str, dict[str, str]] = {}
    if manifest_path:
        path = Path(manifest_path)
        if path.exists():
            with open(path, newline="", encoding="utf-8") as fh:
                manifest = {
                    row["ligand_id"]: row for row in csv.DictReader(fh) if row.get("ligand_id")
                }

    rows = collate(records, manifest, group_stereoisomers=group_stereoisomers)

    with open(run.results_file, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=RESULT_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for rank, row in enumerate(rows, start=1):
            writer.writerow(row.to_row(rank))

    with open(run.all_modes_file, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=ALL_MODES_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for record in records:
            for mode in record.modes:
                writer.writerow(
                    {
                        "Ligand": record.ligand_id,
                        "Mode": mode.mode,
                        "Affinity": mode.affinity,
                        "RMSD_LB": f"{mode.rmsd_lb:.3f}",
                        "RMSD_UB": f"{mode.rmsd_ub:.3f}",
                        "RandomSeed": record.seed if record.seed is not None else "",
                    }
                )

    return run.results_file, run.all_modes_file, len(rows)
