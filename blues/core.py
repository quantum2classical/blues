from rdkit import Chem

from .bondorders import convert_bond_orders
from .buildmol import build_molecule, draw_molecule
from .runjanpa import run_janpa


class MoldenConverter:
    """Converts Molden files into chemical identifier strings."""

    def __init__(self, molden_file):
        """Inits MoldenConverter class."""

        self.molden = molden_file
        self.molecule = None
        self.smiles = None
        self.inchi = None
        self.lewis = None
        return

    def convert(self):
        """Runs JANPA and stores the final molecular structure.""" 
        pass
        return

    def tosmiles(self):
        """Returns SMILES string of molecule in Molden file."""
        pass
        return smiles

    def toinchi(self):
        """Returns InChI string of molecule in Molden file."""
        pass
        return inchi

    def tolewis(self):
        """Returns 2D Lewis structure of molecule in Molden file."""
        pass
        return lewis


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
