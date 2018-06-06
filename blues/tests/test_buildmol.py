import numpy as np
from rdkit import Chem

from ..buildmol import add_atoms, build_molecule


def test_add_atoms():
    atom_list = ['O', 'H', 'H']
    result = add_atoms(atom_list)
    assert isinstance(result, Chem.rdchem.RWMol)
    # assert type(result) == Chem.rdchem.RWMol
    return True


def test_build_molecule():
    atom_list = ['O', 'H', 'H']
    bond_matrix1 = np.array([[2, 0, 'a'], [1, 1, 0], [0, 2, 1]])
    bond_matrix2 = np.array([[1, 1, 0], [0, 2, 1]])
    bond_matrix3 = np.array([[1, 0, 1.2], [0, 2, 0], [0, 0, 1]])
    bond_matrix4 = np.array([[1, 0, 0], [0, 1, 0], [0, 2, 0]])
    # try string value
    try:
        build_molecule(atom_list, bond_matrix1)
    except Exception as err:
        print(err)
    else:
        raise Exception('did not catch string input in bonding matrix')
    # try non-square matrix
    try:
        build_molecule(atom_list, bond_matrix2)
    except Exception as err:
        print(err)
    else:
        raise Exception('did not catch non-square matrix')
    # try float value
    try:
        build_molecule(atom_list, bond_matrix3)
    except Exception as err:
        print(err)
    else:
        raise Exception('did not catch float value not close to an integer')
    # try 0 value diagonal
    try:
        build_molecule(atom_list, bond_matrix4)
    except Exception as err:
        print(err)
    else:
        for i in range(bond_matrix4.shape[0]):
            print(bond_matrix4[i, i])
    # raise Exception('did not catch diagonal value of 0')

    return True
