import os

from PIL import Image

from ..core import MoldenConverter, molden2smiles, molden2inchi, molden2lewis
from ..utils import file_path


def test_MoldenConverter():
    def test_convert():
        test_molden = file_path('test.molden', 'tests/test_files')
        test = MoldenConverter(test_molden)
        
        test.convert()

        real_atoms = ['C', 'C']
        assert test._molecule.


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
