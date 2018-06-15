from rdkit import Chem

from .bondorders import convert_bond_orders
from .buildmol import build_molecule, draw_molecule
from .runjanpa import run_janpa


class MoldenConverter:
    """Converts Molden files into chemical identifier strings."""

    def __init__(self, molden_file):
        """Inits MoldenConverter class.

        Parameters
        ----------
        molden_file : str
            Path to Molden file."""

        self.molden = molden_file
        self._molecule = None
        self.smiles = None
        self.inchi = None
        self.lewis = None
        return

    def tosmiles(self):
        """Returns SMILES string of molecule in Molden file.

        Returns
        -------
        str
            SMILES representation of molecule."""
        if not self._molecule:
            atoms, raw_bonds = run_janpa(self.molden)
            bonds = convert_bond_orders(raw_bonds)
            self._molecule = build_molecule(atoms, bonds)
            self.smiles = Chem.MolToSmiles(self._molecule)
            # two_d, three_d, inchi, smiles = draw_molecule(atoms, bonds)
        else:
            if self.smiles:
                pass
            else:
                self.smiles = Chem.MolToSmiles(self._molecule)

        return self.smiles

    def toinchi(self):
        """Returns InChI string of molecule in Molden file.

        Returns
        -------
        str
            InChI representation of molecule."""
        if not self._molecule:
            atoms, raw_bonds = run_janpa(self.molden)
            bonds = convert_bond_orders(raw_bonds)
            self._molecule = build_molecule(atoms, bonds)
            self.inchi = Chem.MolToInchi(self._molecule)
            # two_d, three_d, inchi, smiles = draw_molecule(atoms, bonds)
        else:
            if self.inchi:
                pass
            else:
                self.inchi = Chem.MolToSmiles(self._molecule)
        return self.inchi

    def tolewis(self):
        """Returns 2D Lewis structure of molecule in Molden file.

        Returns
        -------
        np.ndarray
            Array containing the image of the 2D Lewis structure."""
        if not self.lewis:
            atoms, raw_bonds = run_janpa(self.molden)
            bonds = convert_bond_orders(raw_bonds)
            self._molecule = build_molecule(atoms, bonds)
            self.lewis, three_d, self.inchi, self.smiles = draw_molecule(atoms,
                                                                         bonds)
        else:
            pass
        return self.lewis


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
