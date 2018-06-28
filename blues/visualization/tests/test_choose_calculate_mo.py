from .. import choose_calculate_mo
from choose_calculate_mo import *


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
        
        

def test_my_mo_list():
    
    try:
        my_mo_list(before_homo=10.0, after_homo=23)
    except Exception:
        pass
    else:
        raise Exception('did not catch float value input')
    
    before = 10
    after = -4
    
    orb_string = my_mo_list(10,23)
    assert orb_string == 'homo-10:lumo+-4', 'function produced wrong string input'
    

        
    
    


