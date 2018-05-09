import numpy as np

from ..runjanpa import run_janpa, get_bond_orders, parse_atom_list,\
    parse_bond_orders


def test_run_janpa():

    return


def test_get_bond_orders():

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
