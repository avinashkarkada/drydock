"""Measure docking throughput before committing to a long run.

A screen of a large library is a decision worth making on evidence. The cost per
ligand varies by more than an order of magnitude with box volume, exhaustiveness
and ligand flexibility, so an estimate carried over from a different target is
worth very little.

The measurement that motivated this module: an MMP9 screen whose search box was
derived from nine active-site residues with 5 A padding ran at a median of 81
seconds and a mean of 130 seconds per ligand, with the worst compound taking 888.
That is roughly six days for a 47,000-record library at exhaustiveness 8 on
twelve cores -- a number worth knowing before starting rather than after.

Measure on an idle machine
--------------------------

A benchmark saturates every core it is given, so anything else running competes
with it directly. During development, benchmark runs that overlapped with test
suites and linters on the same twelve cores reported per-ligand costs inflated by
more than 3x, and the discrepancy was initially mistaken for a real effect of a
code change. Two full runs of the same 100 compounds, measured without competing
load, agreed to within 2%.

If a projection here disagrees sharply with what a real screen goes on to do,
suspect the measurement conditions before suspecting the configuration.

Sampling
--------

Ligands are sampled deterministically across the whole directory rather than
taken from the front. Libraries are frequently ordered by something correlated
with size, so a leading slice is systematically unrepresentative -- in CMNPD the
first few hundred compounds are small sulfur-containing molecules that dock an
order of magnitude faster than the median.
"""

from __future__ import annotations

import random
import statistics
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from drydock.core.screen import LigandJob, _dock_one, iter_ligands
from drydock.engines.base import DockConfig

# Sample size.
#
# Raised from 25 after a measured miss: per-ligand cost has a long right tail
# (in a marine natural-product library, a 10x spread between median and worst),
# and total runtime is a sum, so the tail is exactly what the projection depends
# on. A 20-ligand sample of one library gave a mean of 82 s where the full
# 100-compound run gave 122 s -- a 33% underestimate, entirely because the small
# sample missed the slow compounds.
#
# 50 is a compromise: still minutes rather than hours to measure, but wide enough
# that a projection is not dominated by whether one floppy macrocycle happened to
# be drawn.
DEFAULT_SAMPLE = 50


@dataclass(slots=True)
class BenchmarkResult:
    """Measured throughput for one configuration."""

    label: str
    n_sampled: int
    n_ok: int
    n_failed: int
    times: list[float] = field(default_factory=list)
    wall_s: float = 0.0
    box_volume: float | None = None
    exhaustiveness: int | None = None

    @property
    def median_s(self) -> float:
        return statistics.median(self.times) if self.times else 0.0

    @property
    def mean_s(self) -> float:
        return statistics.fmean(self.times) if self.times else 0.0

    @property
    def max_s(self) -> float:
        return max(self.times, default=0.0)

    @property
    def tail_ratio(self) -> float:
        """Mean divided by median: how much the slow tail costs.

        1.0 means a symmetric distribution. Above roughly 1.3 the projection is
        being driven by a minority of slow ligands, and a small sample may well
        have missed them -- so the projection should be read as a lower bound.
        """
        median = self.median_s
        return (self.mean_s / median) if median else 1.0

    def projected_hours(self, n_ligands: int, n_workers: int) -> float:
        """Wall-clock hours to screen a library of this size.

        Uses the mean rather than the median: total time is a sum, and the tail
        of slow ligands contributes to it in a way the median hides. For the run
        that motivated this module, median and mean differed by more than 2x.
        """
        if not self.times or n_workers < 1:
            return 0.0
        return (self.mean_s * n_ligands) / n_workers / 3600

    def to_dict(self) -> dict[str, Any]:
        return {
            "label": self.label,
            "n_sampled": self.n_sampled,
            "n_ok": self.n_ok,
            "n_failed": self.n_failed,
            "median_s": round(self.median_s, 2),
            "mean_s": round(self.mean_s, 2),
            "max_s": round(self.max_s, 2),
            "wall_s": round(self.wall_s, 2),
            "box_volume": self.box_volume,
            "exhaustiveness": self.exhaustiveness,
        }


def sample_ligands(
    pdbqt_dir: str | Path, n: int = DEFAULT_SAMPLE, seed: int = 0
) -> list[tuple[str, Path]]:
    """Pick a representative sample, deterministically.

    Spread across the directory rather than taken from the front, because
    library ordering is rarely independent of molecular size.
    """
    ligands = list(iter_ligands(pdbqt_dir))
    if not ligands:
        raise ValueError(f"no prepared ligands in {pdbqt_dir}")
    if len(ligands) <= n:
        return ligands
    rng = random.Random(seed)
    return [ligands[i] for i in sorted(rng.sample(range(len(ligands)), n))]


def benchmark(
    pdbqt_dir: str | Path,
    config: DockConfig,
    *,
    label: str = "",
    n: int = DEFAULT_SAMPLE,
    n_workers: int = 0,
    seed: int = 0,
) -> BenchmarkResult:
    """Dock a sample and measure how long it took.

    Runs through the same code path as a real screen, so the numbers include
    engine setup and IPC rather than only the docking itself.
    """
    import multiprocessing as mp
    import os

    n_workers = n_workers or (os.cpu_count() or 1)
    sample = sample_ligands(pdbqt_dir, n=n, seed=seed)

    result = BenchmarkResult(
        label=label or config.engine,
        n_sampled=len(sample),
        n_ok=0,
        n_failed=0,
        box_volume=round(config.box.volume, 1),
        exhaustiveness=config.exhaustiveness,
    )

    jobs = [(LigandJob(lid, str(path), want_poses=False), config) for lid, path in sample]

    started = time.perf_counter()
    ctx = mp.get_context("forkserver")
    pool = ctx.Pool(processes=min(n_workers, len(jobs)))
    try:
        for record, _ in pool.imap_unordered(_dock_one, jobs, chunksize=1):
            result.times.append(record.elapsed_s)
            if record.status == "ok":
                result.n_ok += 1
            else:
                result.n_failed += 1
    finally:
        pool.terminate()
        pool.join()

    result.wall_s = time.perf_counter() - started
    return result


def format_comparison(
    results: list[BenchmarkResult], n_ligands: int, n_workers: int
) -> str:
    """Render benchmark results as a table, cheapest first."""
    if not results:
        return "no results"

    lines = [
        f"{'configuration':<28} {'volume':>9} {'exh':>4} "
        f"{'median':>8} {'mean':>8} {'max':>8} {'projected':>11}",
        "-" * 82,
    ]
    for r in sorted(results, key=lambda r: r.mean_s):
        volume = f"{r.box_volume:,.0f}" if r.box_volume else "-"
        hours = r.projected_hours(n_ligands, n_workers)
        projected = f"{hours:.1f} h" if hours < 48 else f"{hours / 24:.1f} d"
        lines.append(
            f"{r.label:<28} {volume:>9} {r.exhaustiveness or '-':>4} "
            f"{r.median_s:>7.1f}s {r.mean_s:>7.1f}s {r.max_s:>7.1f}s {projected:>11}"
        )

    lines.append("")
    lines.append(
        f"projection assumes {n_ligands:,} ligands across {n_workers} workers, "
        "scaled from the mean"
    )

    heavy = [r for r in results if r.tail_ratio > 1.3]
    if heavy:
        worst = max(r.tail_ratio for r in heavy)
        lines.append(
            f"note: cost is skewed by slow ligands (mean/median up to {worst:.1f}x). "
            "A sample can miss those, so read these projections as lower bounds."
        )
    return "\n".join(lines)
