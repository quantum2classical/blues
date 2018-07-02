from .. import choose_calculate_mo
from choose_calculate_mo import get_atom_positions
from choose_calculate_mo import *
from ...utils import file_path


def test_read_molden():
    
    try:
        read_molden(filename='nw.cube', itype='molden')
    except Exception:
        pass
    else:
        raise Exception('did not ensure filetype matches itype')
    
    try:
        read_molden(filename=2)
    except Exception:
        pass
    else:
        raise Exception('did not catch wrong input type')
    
    return True
        
        

def test_my_mo_list():
    
    try:
        my_mo_list(before_homo=10.0, after_homo=23)
    except Exception:
        pass
    else:
        raise Exception('did not catch float value input')
    
    before = 10
    after = -4
    
    orb_string = my_mo_list(before, after)
    assert orb_string == 'homo-10:lumo+-4', 'function produced wrong string input'
    
    before = 10
    after = -20
    
    try:
        orb_string = my_mo_list(before, after)
    except Exception:
        pass
    else:
        raise Exception('did not catch case of after<0 and abs(after)>abs(before)')
    
    return True
    
def test_get_atom_positions():
    # get path to test file
    pth = file_path('test_nw.molden', 'visualization/tests/test_files')
    
    xyz_c, xyz_o, xyz_n = get_atom_position(pth)
    
    assert len(xyz_c[0]) == 21, 'failed to parse correct number of C atoms'
    assert xyz_c[0][2] == 6.403137996219281, 'error in parsing test_file'
    
    return True
    
def test_calculate_densities():
    # get path to test file
    pth = file_path('test_nw.molden', 'visualization/tests/test_files')
    
    x, y, z, mo_list, mo_info = calculate_densities(pth, before_homo=1, 
                                after_homo=0)
    assert type(mo_info) == dict
    
    return True
    
    
    
    
    
    



    

    

    

        
    
    


