import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw


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