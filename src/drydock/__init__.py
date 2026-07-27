"""Drydock: reproducible high-throughput virtual screening.

Three stages, kept separate:

``prep-receptor``
    Read a PDB or mmCIF, add hydrogens and charges, and check the result.

``prep-ligands``
    Turn a compound library (SDF/SMILES/MOL2) into docking-ready PDBQT files
    plus a descriptor manifest.

``screen``
    Dock a prepared library against the receptor and emit a ranked CSV.

Zinc pseudo-atom placement for AutoDock4Zn is opt-in and lives in
:mod:`drydock.core.zinc`, since the AD4Zn scoring path cannot run without it.
"""

__version__ = "0.1.0"

__all__ = ["__version__"]
