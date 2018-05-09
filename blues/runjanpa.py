import numpy as np
from subprocess import call


def run_janpa(molden_file):
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

    return (elements, bond_orders)


def parse_atom_list(line):
    """Converts single string to list of elemental symbols.

    Parameters
    ----------
    line : str

    Returns
    -------
    atom_list : list
    """

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

    return bond_orders
