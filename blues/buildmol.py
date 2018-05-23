import os
#from IPython.display import Image
#import numpy as np
#from rdkit import Chem
#from rdkit.Chem import AllChem, Draw

def order_atoms(filename='atom_list.txt'):
    """list out the txt file containing all atoms in the molecule, order the same as JANPA output and eliminate whitespace"""
    with open(filename) as f:
        atoms = f.readlines()
        atoms = [atom.strip() for atom in atoms]

    return atoms


def add_atoms(atoms):
    """add atoms and correct bonds to build molecule"""
    # create molecule under RWMol class, easily editable
    mol = Chem.RWMol()
    # call function to order atoms
    order_atoms()
    # add all atoms to molecule
    for atom in atoms:
        mol.AddAtom(Chem.Atom(atom))
    return mol

def build_molecule(bonding_matrix):
    """Takes connection matrix and builds the molecule with RDKit"""
    # add bonds from connection matrix built from bondorders.py
    bonds = bonding_matrix
    # ensure matrix is square
    assert bonding_matrix.shape[0] == bonding_matrix.shape[1], 'matrix must be square!'
    # ensure diagonals don't equal 0
    for i in range(bonding_matrix.shape[0]):
        assert bonding_matrix[i,i] != 0, 'diagonals cannot equal 0'
    for i, j in zip(*np.triu_indices_from(bonds, 1)):
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


def draw_molecule(mol):
    """represent built molecule as 2 and 3 dimensional images"""
    # use RDKit internal check to verify molecule is physical
    Chem.SanitizeMol(mol)
    mol_noH = Chem.RemoveHs(mol)
    # 2d representation of molecule
    two_d = Draw.MolToImage(mol_noH)
    # 3d picture
    three_d = Image('mol.png', width=400)
    inchi = "InChI: {}".format(Chem.MolToInchi(mol_noH))
    smiles = "SMILES: {}".format(Chem.MolToSmiles(mol_noH))

    return two_d, three_d, inchi, smiles