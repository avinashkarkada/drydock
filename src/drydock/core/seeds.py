"""Deterministic per-ligand seeds.

Both ligand preparation (ETKDG embedding) and docking (Vina's Monte Carlo search)
are stochastic. Reproducing a run therefore means reproducing the seed each
ligand was given.

Deriving seeds by counting, ``global_seed + i``, would tie them to the order
work happened to be scheduled in, and a parallel run does not have a stable
order. Instead each seed is a hash of the global seed and the ligand's
identifier, so a ligand gets the same seed regardless of how many workers ran,
which ligands failed, or whether the run was resumed halfway through.

Sequential seeds also correlate: pseudo-random generators started at N and N+1
can produce related streams. Hashing decorrelates them for free.
"""

from __future__ import annotations

import hashlib

# AutoDock Vina takes a signed 32-bit seed, which is the tightest constraint of
# anything consuming these, so everything is generated in that range.
_SEED_MAX = 2**31 - 1


def ligand_seed(global_seed: int, ligand_id: str) -> int:
    """Derive a stable seed for one ligand.

    The same ``(global_seed, ligand_id)`` pair always yields the same value,
    across processes, machines and Python versions, ``hash()`` is not used
    because that is randomised per process.

    Args:
        global_seed: The run's seed, recorded in provenance.
        ligand_id: Unique ligand identifier.

    Returns:
        A seed in ``[0, 2**31 - 1]``.
    """
    digest = hashlib.sha256(f"{global_seed}:{ligand_id}".encode()).digest()
    return int.from_bytes(digest[:4], "big") % _SEED_MAX


def conformer_seed(global_seed: int, ligand_id: str, conformer: int) -> int:
    """Derive a stable seed for one conformer of one ligand."""
    return ligand_seed(global_seed, f"{ligand_id}#conf{conformer}")
