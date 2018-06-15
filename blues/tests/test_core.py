import os

from PIL import Image

from ..core import MoldenConverter, molden2smiles, molden2inchi, molden2lewis
from ..utils import file_path


# def test_MoldenConverter():
#    # Define all variables and paths to be used in method tests.
#    test_molden = file_path('test.molden', 'tests/test_files')
#    real_atoms = ['C', 'C', 'C', 'C', 'C', 'C',
#                  'C', 'C', 'C', 'C', 'C', 'C']
#    real_smiles = 'C=CCc1cccc(C#CC)c1'
#    real_inchi =\
#        'InChI=1S/C12H12/c1-3-6-11-8-5-9-12(10-11)7-4-2/h3,5,8-10H,1,6H2,2H3'
#    real_lewis_path = file_path('test_mol.png', 'tests/test_files')
#    real_lewis = Image.open(real_lewis_path)


def test_tosmiles():
    test_molden = file_path('test.molden', 'tests/test_files')
    real_atoms = ['C', 'C', 'C', 'C', 'C', 'C',
                  'C', 'C', 'C', 'C', 'C', 'C']
    real_smiles = 'C=CCc1cccc(C#CC)c1'
    print('testing smiles')
    test = MoldenConverter(test_molden)
    test_smiles = test.tosmiles()

    assert test._molecule.GetNumAtoms() == 12
    assert [atom.GetSymbol() for atom
            in test._molecule.GetAtoms()] == real_atoms
    assert test_smiles == real_smiles
    return


def test_toinchi():
    test_molden = file_path('test.molden', 'tests/test_files')
    real_atoms = ['C', 'C', 'C', 'C', 'C', 'C',
                  'C', 'C', 'C', 'C', 'C', 'C']
    real_inchi =\
        'InChI=1S/C12H12/c1-3-6-11-8-5-9-12(10-11)7-4-2/h3,5,8-10H,1,6H2,2H3'
    print('testing inchi')
    test = MoldenConverter(test_molden)
    test_inchi = test.toinchi()

    assert test._molecule.GetNumAtoms() == 12
    assert [atom.GetSymbol() for atom
            in test._molecule.GetAtoms()] == real_atoms
    assert test_inchi == real_inchi
    return


def test_tolewis():
    test_molden = file_path('test.molden', 'tests/test_files')
    real_atoms = ['C', 'C', 'C', 'C', 'C', 'C',
                  'C', 'C', 'C', 'C', 'C', 'C']
    real_lewis_path = file_path('test_mol.png', 'tests/test_files')
    real_lewis = Image.open(real_lewis_path)
    print('testing lewis')
    test = MoldenConverter(test_molden)
    test_lewis = test.tolewis()
    test_lewis.save('test_lewis.png')
    test_lewis = Image.open('test_lewis.png')

    assert test._molecule.GetNumAtoms() == 12
    assert [atom.GetSymbol() for atom
            in test._molecule.GetAtoms()] == real_atoms
    try:
        assert test_lewis == real_lewis, "Images do not match."
    except(AssertionError):
        os.remove('test_lewis.png')
        raise AssertionError("Images do not match.")
    else:
        os.remove('test_lewis.png')

    return


def test_molden2smiles():
    real_smiles = 'C=CCc1cccc(C#CC)c1'
    test_molden = file_path('test.molden', 'tests/test_files')

    test_smiles = molden2smiles(test_molden)

    assert test_smiles == real_smiles

    return


def test_molden2inchi():
    real_inchi =\
        'InChI=1S/C12H12/c1-3-6-11-8-5-9-12(10-11)7-4-2/h3,5,8-10H,1,6H2,2H3'
    test_molden = file_path('test.molden', 'tests/test_files')

    test_inchi = molden2inchi(test_molden)

    assert test_inchi == real_inchi

    return


def test_molden2lewis():
    real_lewis_path = file_path('test_mol.png', 'tests/test_files')
    real_lewis = Image.open(real_lewis_path)
    test_molden = file_path('test.molden', 'tests/test_files')

    test_lewis = molden2lewis(test_molden)
    test_lewis.save('test_lewis.png')
    test_lewis = Image.open('test_lewis.png')

    try:
        assert test_lewis == real_lewis, "Images do not match."
    except(AssertionError):
        os.remove('test_lewis.png')
        raise AssertionError("Images do not match.")
    else:
        os.remove('test_lewis.png')

    return
