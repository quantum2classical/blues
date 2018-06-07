import re
import os

import numpy as np
import subprocess
from subprocess import call


from .utils import file_path


def run_janpa(molden_file, wiberg_file='wiberg.dat', delete_wiberg_file=False):
    """Runs JANPA with specified input file.

    Parameters
    ----------
    molden_file : str
        Name of MOLDEN file automatically generated from NWChem run.

    Returns
    -------
    elements : list
        List of atoms in system in the same order as input file.
    bond_orders : np.ndarray
        Wiberg-Mayer bond orders calculated in JANPA.
    """
    if not isinstance(molden_file, str):
        raise TypeError("Argument must be a string.")
    elif not os.path.exists(molden_file):
        raise OSError("Input file not found.")
    else:
        pass

    janpa_path = file_path('janpa.jar', 'janpa')

    retcode = call(['java', '-jar', janpa_path, '-i', molden_file,
                    '-ignoreFock', '-WibergBondOrders_File', wiberg_file],
                   stdout=subprocess.DEVNULL)

    if retcode != 0:
        raise Exception("JANPA program did not exit with code 0.")
    else:
        pass

    elements, bond_orders = get_bond_orders(wiberg_file)

    if delete_wiberg_file:
        os.remove(wiberg_file)
    else:
        pass

    return (elements, bond_orders)


def get_bond_orders(wiberg_file):
    """Extracts the Wiberg-Mayer bond orders from generated file.

    Parameters
    ----------
    wiberg_file : str
        Name of file containing Wiberg-Mayer bond orders generated
        from JANPA execution.

    Returns
    -------
    elements : list
        List of atoms in system in the same order as input file.
    bond_orders : np.ndarray
        Wiberg-Mayer bond orders calculated in JANPA.
    """
    if not isinstance(wiberg_file, str):
        raise TypeError("File containing Wiberg bond "
                        "orders must be a string.")
    else:
        pass

    with open(wiberg_file, 'r') as f:
        f.readline()  # skip header explaining Wiberg-Mayer bond indices
        num_atoms = f.readline()  # read in shape of matrix
        atom_list = f.readline()  # get atom list
        f.readline()  # skip blank line
        matrix = f.readlines()  # get Wiberg-Mayer bond orders

    matrix_shape = tuple(int(num) for num in num_atoms.split())

    elements = parse_atom_list(atom_list)
    bond_orders = parse_bond_orders(matrix)

    assert bond_orders.shape == matrix_shape,\
        "Incorrect bond order matrix shape."

    return (elements, bond_orders)


def parse_atom_list(raw_atom_list):
    """Converts single string to list of elemental symbols.

    Parameters
    ----------
    raw_atom_list : str

    Returns
    -------
    atom_list : list
    """

    split_line = raw_atom_list.split()
    for group in split_line:
        if not group.istitle():
            raise ValueError("Input must be atomic symbols "
                             "with first letter capitalized.")
        else:
            pass

    element = re.compile('[A-Z][a-z]?')

    # Splits each line into separate fields and pulls the element symbol.
    atom_list = [element.match(atom).group() for atom in raw_atom_list.split()]

    return atom_list


def parse_bond_orders(raw_bond_orders):
    """Converts a list of string to a matrix of bond orders.

    Parameters
    ----------
    raw_bond_orders : list of str

    Returns
    -------
    bond_orders : np.ndarray
    """
    # Splits each line into separate fields and converts to float,
    # ignores last column that contains atom identifier.
    bond_orders = np.array([[float(num) for num in line.split()[:-1]]
                            for line in raw_bond_orders])
    i, j = np.tril_indices(len(bond_orders))
    bond_orders[i, j] = 0
    return bond_orders
