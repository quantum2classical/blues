from rdkit import Chem

from .bondorders import convert_bond_orders
from .buildmol import build_molecule, draw_molecule
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


def molden2inchi(molden_file):
    """Uses molden file as input to JANPA and generates an InChI string.

    Parameters
    ----------
    molden_file : str

    Returns
    -------
    inchi : str
    """
    atoms, raw_bonds = run_janpa(molden_file)
    bonds = convert_bond_orders(raw_bonds)
    mol = build_molecule(atoms, bonds)
    inchi = Chem.MolToInchi(mol)

    return inchi


def molden2lewis(molden_file, filename=None, save=False):
    """Uses molden file as input to JANPA and generates a 2D Lewis structure.

    Parameters
    ----------
    molden_file : str

    Returns
    -------
    lewis : image
    """
    atoms, raw_bonds = run_janpa(molden_file)
    bonds = convert_bond_orders(raw_bonds)
    lewis, three_d, inchi, smiles = draw_molecule(atoms, bonds)

    return lewis
