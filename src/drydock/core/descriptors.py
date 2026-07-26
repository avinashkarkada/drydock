"""Molecular descriptors for the ligand manifest.

These exist to make a hit list interpretable. A screen that returns
``CMNPD31204, -11.8`` is not actionable; the same row with a molecular weight, a
logP and a ligand efficiency tells you whether the score is worth chasing or is
just a large greasy molecule scoring well because it is large and greasy.

Descriptors are computed on the **input** molecule -- sanitised, before
protonation -- because that is the form people expect library statistics to
describe. Properties of the prepared, protonated species that actually gets
docked are recorded separately.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

# Ligand efficiency divides binding energy by heavy-atom count. It is the standard
# correction for the fact that docking scores grow with molecular size, which
# matters especially for natural-product libraries where the mass range is wide.
# Reported in kcal/mol per heavy atom.


@dataclass(frozen=True, slots=True)
class Descriptors2D:
    """Per-compound descriptors, one row of the ligand manifest."""

    smiles: str
    formula: str
    mw: float
    heavy_atoms: int
    rot_bonds: int
    clogp: float
    tpsa: float
    hbd: int
    hba: int
    formal_charge: int
    rings: int
    aromatic_rings: int
    fraction_csp3: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def compute(mol: Chem.Mol) -> Descriptors2D:
    """Compute descriptors for a sanitised RDKit molecule.

    Hydrogens are stripped first so counts describe the compound rather than the
    particular hydrogen treatment of whatever file it arrived in.
    """
    bare = Chem.RemoveHs(mol) if mol.GetNumAtoms() else mol

    return Descriptors2D(
        smiles=Chem.MolToSmiles(bare),
        formula=rdMolDescriptors.CalcMolFormula(bare),
        mw=round(Descriptors.MolWt(bare), 3),
        heavy_atoms=bare.GetNumHeavyAtoms(),
        rot_bonds=Descriptors.NumRotatableBonds(bare),
        clogp=round(Crippen.MolLogP(bare), 3),
        tpsa=round(Descriptors.TPSA(bare), 2),
        hbd=Descriptors.NumHDonors(bare),
        hba=Descriptors.NumHAcceptors(bare),
        formal_charge=Chem.GetFormalCharge(bare),
        rings=rdMolDescriptors.CalcNumRings(bare),
        aromatic_rings=rdMolDescriptors.CalcNumAromaticRings(bare),
        fraction_csp3=round(rdMolDescriptors.CalcFractionCSP3(bare), 4),
    )


def ligand_efficiency(affinity: float | None, heavy_atoms: int | None) -> float | None:
    """Binding energy per heavy atom.

    Returns None rather than raising on missing inputs, because a failed ligand
    legitimately has no affinity and the manifest may not cover every compound in
    a resumed run.
    """
    if affinity is None or not heavy_atoms:
        return None
    return round(affinity / heavy_atoms, 4)
