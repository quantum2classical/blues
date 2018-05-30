from ..core import molden2smiles
from ..utils import file_path


def test_molden2smiles():

    real_smiles = 'C=CCc1cccc(C#CC)c1'
    test_molden = file_path('test.molden', 'tests/test_files')

    test_smiles = molden2smiles(test_molden)

    assert test_smiles == real_smiles
    
    return
