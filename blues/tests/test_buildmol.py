from ..buildmol import *

def test_order_atoms():
    """   """
    order_atoms('atom_list.txt')
    assert sum(line.isspace() for line in atoms) == 0, 'order_atoms function failed to eliminate whitespace'
    actual_length = len(atoms)
    assert len(atoms) == 47, 'molecule should have 47 atoms,' + str(actual_length) + 'was counted'
    return atoms


def test_add_atoms():
    add_atoms(atoms)