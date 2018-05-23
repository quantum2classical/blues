import numpy as np

from ..runjanpa import run_janpa, get_bond_orders, parse_atom_list,\
    parse_bond_orders
from ..utils import file_path


def test_run_janpa():
    try:
        run_janpa('fake.molden')
    except(OSError):
        pass
    else:
        raise Exception("Function executed with nonexistent input file.")

    try:
        run_janpa(1234)
    except(TypeError):
        pass
    else:
        raise Exception("Function executed with numerical argument.")

    test_molden_file = file_path('test_files/test.molden', 'tests')

    atoms, bonds = run_janpa(test_molden_file, delete_wiberg_file=True)

    real_atoms = ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'C', 'H', 'H', 'H',
                  'C', 'C', 'H', 'H', 'H', 'C', 'H', 'C', 'H', 'C', 'H',
                  'H', 'H']

    aromatic_bonds = (np.array([0, 0, 1, 2, 3, 4]),
                      np.array([1, 2, 3, 5, 4, 5]))

    double_bond = (np.array([18]), np.array([20]))
    triple_bond = (np.array([7]), np.array([11]))

    assert atoms == real_atoms, "Incorrect atoms extracted."
    aromatic_id = np.where(np.isclose(bonds, 1.5, atol=0.20))
    assert np.all(aromatic_id[0] == aromatic_bonds[0])
    assert np.all(aromatic_id[1] == aromatic_bonds[1])
    assert np.where(np.isclose(bonds, 2.0, atol=0.30)) == double_bond
    assert np.where(np.isclose(bonds, 3.0, atol=0.30)) == triple_bond

    return


def test_get_bond_orders():

    try:
        get_bond_orders('fake.txt')
    except(OSError):
        pass
    else:
        raise Exception("Function executed with nonexistent input file.")

    try:
        get_bond_orders(1234)
    except(TypeError):
        pass
    else:
        raise Exception("Function executed with numerical argument.")

    test_wiberg_file = file_path("test_files/wiberg.txt", 'tests')
    atoms, bonds = get_bond_orders(test_wiberg_file)

    real_atoms = ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'C', 'H', 'H', 'H',
                  'C', 'C', 'H', 'H', 'H', 'C', 'H', 'C', 'H', 'C', 'H',
                  'H', 'H']

    aromatic_bonds = (np.array([0, 0, 1, 2, 3, 4]),
                      np.array([1, 2, 3, 5, 4, 5]))

    double_bond = (np.array([18]), np.array([20]))
    triple_bond = (np.array([7]), np.array([11]))

    assert atoms == real_atoms, "Incorrect atoms extracted."
    aromatic_id = np.where(np.isclose(bonds, 1.5, atol=0.20))
    assert np.all(aromatic_id[0] == aromatic_bonds[0])
    assert np.all(aromatic_id[1] == aromatic_bonds[1])
    assert np.where(np.isclose(bonds, 2.0, atol=0.30)) == double_bond
    assert np.where(np.isclose(bonds, 3.0, atol=0.30)) == triple_bond

    return


def test_parse_atom_list():

    test_line1 = "c1\tc2\tc3\tc4\tc5\tc6\th7\th8\th9\th10\th11\th12\t\n"
    try:
        parse_atom_list(test_line1)
    except(ValueError):
        pass
    else:
        raise Exception("Bad input allowed.")

    test_line2 = "C1\tC2\tC3\tC4\tC5\tC6\tH7\tH8\tH9\tH10\tH11\tH12\t\n"
    test_out = parse_atom_list(test_line2)

    assert test_out == ['C'] * 6 + ['H'] * 6, "Atom list is incorrect."

    return


def test_parse_bond_orders():

    test1 = ['0.1\t0.2\t0.3\t\tH1\n',
             '0.4\t0.5\t0.6\t\tH2\n',
             '0.7\t0.8\t0.9\t\tO3\n']

    test_out = parse_bond_orders(test1)

    assert isinstance(test_out, np.ndarray), "Output is not an ndarray."
    assert test_out.dtype == np.float64, "dtype is not float64."
    assert test_out.shape[0] == test_out.shape[1], "Output is not square."

    return
