import os
from IPython.display import Image
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw



def convert_bond_orders(bond_orders):
## Round bond orders to their nearest bond type
    for i in range(bonds.shape[0]):
        for j in range(bonds.shape[1]):
            if np.isclose(bonds[i, j], 1.5, atol=0.15):
                bonds[i, j] = 1.5
            else:
                bonds[i, j] = int(round(bonds[i, j]))
    bonding_matrix = bounds
    return bonding_matrix