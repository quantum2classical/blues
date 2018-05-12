import os
from IPython.display import Image
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

def order_atoms(filename='atom_list.txt')
    """read out the txt file containing all atoms in the molecule, order the same as JANPA output and eliminate whitespace"""
    with open(filename) as f:
        atoms = f.readlines()
        atoms = [atom.strip() for atom in atoms]

mol = Chem.RWMol()

def
for atom in atoms:
    mol.AddAtom(Chem.Atom(atom))

for i, j in zip(*np.triu_indices_from(bonds, 1)):
    if np.isclose(bonds[i, j], 1):
        mol.AddBond(int(i), int(j), Chem.BondType.SINGLE)
    elif np.isclose(bonds[i, j], 1.5):
        mol.AddBond(int(i), int(j), Chem.BondType.AROMATIC)
    elif np.isclose(bonds[i, j], 2):
        mol.AddBond(int(i), int(j), Chem.BondType.DOUBLE)
    elif np.isclose(bonds[i, j], 3):
        mol.AddBond(int(i), int(j), Chem.BondType.TRIPLE)
    else:
        pass

for at in mol.GetAtoms():
    at.SetNoImplicit(True)

Chem.SanitizeMol(mol)

# 2d representation of molecule
mol_noH = Chem.RemoveHs(mol)
Draw.MolToImage(mol_noH)