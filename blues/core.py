from rdkit import Chem

from .bondorders import convert_bond_orders
from .buildmol import build_molecule
from .runjanpa import run_janpa


def molden2smiles(molden_file):
    """
    Uses molden file as input to JANPA and generates a smiles string.

    Parameters
    ----------
    molden_file : str

    Returns
    -------
    smiles : str
    """

    atoms, raw_bonds = run_janpa(molden_file)
    bonds = convert_bond_orders(raw_bonds)
    mol = build_molecule(atoms, bonds)
    smiles = Chem.MolToSmiles(mol)
    # two_d, three_d, inchi, smiles = draw_molecule(atoms, bonds)

    return smiles
