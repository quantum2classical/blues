import os
from IPython.display import Image
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw

def add_atoms(atom_file):
    """
    Adds atoms and correct bonds to build molecule using an RDKit class

    Parameters
    ----------
    atom_file: string of file path

    Returns: RWMol object

    """
    # create molecule under RWMol class
    mol = Chem.RWMol()
    # add all atoms to molecule
    for atom in atom_file:
        mol.AddAtom(Chem.Atom(atom))
    return mol



def build_molecule(atom_list, bo_matrix):
    """
    Takes .txt file, builds molecule and connection matrix specifying bonds, and builds the molecule with RDKit

    Parameters
    ----------
    atom_file: string of file path

    Returns:

    RWMol object with structure and bonding based on JANPA's Wiberg-Mayer bond index matrix

    """
    # ensure matrix is square
    assert bo_matrix.shape[0] == bo_matrix.shape[1], 'matrix must be square!'
    # ensure diagonals don't equal 0
    for i in range(bo_matrix.shape[0]):
        assert bo_matrix[i,i] != 0, 'diagonals cannot equal 0'
    for i, j in zip(*np.triu_indices_from(bo_matrix, 1)):
        if np.isclose(bonds[i, j], 1):
            mol.AddBond(int(i), int(j), Chem.BondType.SINGLE)
        elif np.isclose(bonds[i, j], 1.5):
            mol.AddBond(int(i), int(j), Chem.BondType.AROMATIC)
        elif np.isclose(bonds[i, j], 2):
            mol.AddBond(int(i), int(j), Chem.BondType.DOUBLE)
        elif np.isclose(bonds[i, j], 3):
            mol.AddBond(int(i), int(j), Chem.BondType.TRIPLE)
        elif np.isclose(bonds[i, j], 0):
            pass
        else:
            raise Exception('indeterminate or nonphysical bond type!')
    # Do not show hydrogens
    mol = add_atoms(atoms)
    for at in mol.GetAtoms():
        at.SetNoImplicit(True)
    return mol


def draw_molecule(atom_list):
    """
    Represents built molecule as 2 and 3 dimensional images and output INCHI and SMILES. Molecule is built from
    txt file containing ordered elements, bonds are added from Wiberg-Mayer matrix contained in JANPA output

    Parameters
    ----------
    atom_list: string containing file path

    Returns:
    ----------
    2d and 3d images as well as INCHI and SMILES strings for the molecule
    """
    # build molecule
    mol = build_molecule(atom_list)
    # 2d representation of molecule
    two_d = Draw.MolToImage(mol)
    # 3d picture
    three_d = Image('mol.png', width=400)
    inchi = "InChI: {}".format(Chem.MolToInchi(mol))
    smiles = "SMILES: {}".format(Chem.MolToSmiles(mol))

    return two_d, three_d, inchi, smiles