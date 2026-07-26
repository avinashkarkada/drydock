"""Drydock: reproducible high-throughput virtual screening.

Two stages, deliberately decoupled:

``prep-ligands``
    Turn a compound library (SDF/SMILES/MOL2) into docking-ready PDBQT files
    plus a descriptor manifest.

``screen``
    Dock a prepared library against a prepared receptor and emit a ranked CSV.

Receptor preparation is out of scope: point Drydock at a receptor PDBQT you
prepared yourself. The one exception is zinc pseudo-atom placement, which is
opt-in and lives in :mod:`drydock.core.zinc`, because the AutoDock4Zn scoring
path is unusable without it.
"""

__version__ = "0.1.0"

__all__ = ["__version__"]
