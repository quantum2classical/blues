from ..bondorders import *

def test_convert_bond_orders():
    bond_orders = get_bond_orders('output.janpa')##should be modified
    bonding_matrix = convert_bond_orders(bond_orders)
    assert bonding_matrix.shape[0] == bonding_matrix.shape[1], "the bond orders matrix should be square"
    assert bonding_matrix[2][2] > 0, "the diagonal of the bond orders matrix should be positive"
    return