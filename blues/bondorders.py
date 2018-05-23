import os
from IPython.display import Image
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

def tail(f, lines=1, _buffer=4098):
    """Tail a file and get X lines from the end"""
    # place holder for the lines found
    lines_found = []

    # block counter will be multiplied by buffer
    # to get the block size from the end
    block_counter = -1

    # loop until we find X lines
    while len(lines_found) < lines:
        try:
            f.seek(block_counter * _buffer, os.SEEK_END)
        except IOError:  # either file is too small, or too many lines requested
            f.seek(0)
            lines_found = f.readlines()
            break

        lines_found = f.readlines()

        # decrement the block counter to get the
        # next X bytes
        block_counter -= 1

    return lines_found[-lines:]


def get_bond_orders(filename):
    """Fetches the bond indices from the end of the JANPA output."""
    f = open(filename)
    lines = tail(f, lines=56)
    f.close()

    split_lines = [line.split() for line in lines[:-9]]

    for line in split_lines:
        del line[0:2]
        line[0] = line[0][:-1]

    n_atoms = len(split_lines)
    bond_orders = np.zeros((n_atoms, n_atoms))

    indices = np.triu_indices(47)

    bond_orders[indices] = [float(item) for sublist in split_lines for item in sublist]
    return bond_orders

def convert_bond_orders(bond_orders):
    """Round bond orders to their nearest bond type"""
    for i in range(bond_orders.shape[0]):
        for j in range(bond_orders.shape[1]):
            if np.isclose(bond_orders[i, j], 1.5, atol=0.15):
                bond_orders[i, j] = 1.5
            else:
                bond_orders[i, j] = int(round(bond_orders[i, j]))
    bonding_matrix = bond_orders
    return bonding_matrix