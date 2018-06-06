import numpy as np

from ..bondorders import convert_bond_orders


def test_convert_bond_orders():
    bond_orders = np.random.random((5, 5)) * 0.01 + 1.5
    assert convert_bond_orders(bond_orders).shape[0] ==\
        convert_bond_orders(bond_orders).shape[1],\
        "the bond orders matrix should be square"
    assert convert_bond_orders(bond_orders)[2][2] == 1.5,\
        "all values about 1.5 should be reduced into exact 1.5"
    bond_orders = np.random.random((5, 5)) * 0.1 + 10
    assert convert_bond_orders(bond_orders)[2][2] == 10,\
        "other values aren't close to 1.5 should be reduced into integer"
    return
