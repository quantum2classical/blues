import numpy as np
from . import *
import pandas as pd
import sys
sys.path.append('C:\orbkit-cython')
from orbkit import read, grid, display, atomic_populations, options, main, extras, output
from orbkit.output import pdb_creator ,xyz_creator
import orbkit as ok
# this submodule should output:
# 1. strings of the number of MOs
# 2. a list of which MOs will be calculated
# 3. calculated density plots for all MOs
# 4. 3d grid with x, y, and z values



def read_molden(filename, all_mo=True):
"""
string of molden file, calculate all mo's

"""
    # read molden file into qc class
    read_file = read.main_read(filename, itype='molden', all_mo=True)
    return read_file


def enumerate_mo(filename, all_mo=True):
"""
input: molden file name
returns: molden file read by orbkit,
        string telling how many MOs there are, which is HOMO
"""
    read_file = read_molden(filename)
    #find out how many MOs and where HOMO is by finding when MO occupancy = 0
    for x in range(len(qc.mo_spec)):
        if qc.mo_spec[x]['occ_num'] != 0.0:
            continue
        else:
            where_is_homo = print('HOMO is MO number ' + str(x - 1) + \
              ' out of ' + str(len(qc.mo_spec)))
            break
            
    print(where_is_homo)
    
    return read_file


def my_mo_list(filename, before_homo, after_homo, all_mo=True):
"""
input: molden filename, how many MOs to calculate before HOMO, 
how many after HOMO

returns:
"""
    # 
    read_file = enumerate_mo(filename)
    
    # make string to input into orbkit's calculation function
    orbital_string_input = 'homo-'+str(before_homo)+':lumo+'+str(after_homo)
    
    # make list of orbitals user wishes to have calculated
    orbital_list = []
    orb_list = ['HOMO+'+str(x) for x in range(-before_homo, after_homo+1)]

    for x in range(len(orb_list)):
        if '-' in orb_list[x]:
            orbital_list.append(orb_list[x].replace('+', ''))
        elif '0' in orb_list[x]:
            orbital_list.append(orb_list[x].replace('+0',''))
        elif orb_list[x]=='HOMO+1':
            orbital_list.append('LUMO')
        else:
            orbital_list.append(orb_list[x])
       
    print(orbital_list)
     
    return read_file, orbital_string_input



def get_atom_positions(filename, output_name='charges', comments=''):
"""
parameters
----------
filename: 
    
returns
--------
PDB file with atomic positions and  

"""
    # read molden file, create qc class
    read_file = read_molden(filename, itype='molden', all_mo=True, 
    output_name='charges', comments='') 
    
    # create dictionary of mulliken, lowdin pop analysis
    pop = atomic_populations.mulliken(read_file)
    pop = atomic_populations.lowdin(read_file)
    # output pdb file to get atomic positions and partial charges
    pdb_creator(read_file.geo_info, read_file.geo_spec, 
    output_name=output_name, charges=pop['charge'],comments=comments) 
     
    # read file
    atoms = pd.read_csv(output_name+'.pdb',header=None, delim_whitespace=True,
    names=['what','number','atom','x','y','z','charge'])
    
    # parse file to obtain positions
    atoms=data.loc[3:len(data)-3]
    atoms=data.iloc[:,2:6]
    # remove hydrogens
    atoms = data[data.atom != 'H']
    # create individual dataframes for O, N atoms
    atoms_o = data[data.atom=='O']
    atoms_n = data[data.atom =='N']
    # remove O,N from original dataframe
    atoms = data[data.atom != 'O']
    atoms = data[data.atom !='N'] 
    
    # convert all to arrays
    xc = np.array(atoms.x)
    yc = np.array(atoms.y)
    zc = np.array(atoms.z)
    xo = np.array(atoms_o.x)
    yo = np.array(atoms_o.y)
    zo = np.array(atoms_o.z)
    xn = np.array(atoms_n.x)
    yn = np.array(atoms_n.y)
    zn = np.array(atoms_n.z)
    
    xyz_c = [xc, yc, zc]
    xyz_o = [xo, yo, zo]
    xyz_n = [xn, yn, zn]
    
    return xyz_c, xyz_o, xyz_n
    

        
def calculate_densities(filename, before_homo, after_homo, 
                    extend=7.0, output_name='charges', comments='', 
                    density_done=False)
"""
parameters
---------
molden filename: str
before_homo: integer of first MO you want to calcualate
after_homo: integer of last MO you want to calculate
extend: float of how much to extend grid so MOs can be visualized. 7.0 default

returns
---------
where_is_homo: string of number of MOs and where HOMO lies
orbital_list: list of MOs that will be calculated
x, y, z: arrays making up grid
mo_list: list of dictionaries containing MO information
mo_info: dict

"""

global density_done
    read_file,  orbital_string_input 
    = my_mo_list(filename, before_homo, after_homo)
    
    # create file with atom positions
    
    # adjust grid to fit molecule + extended units
    # this prevents calc_mo function from automatically setting up grid     
    grid.adjust_to_geo(qc,extend=extend,step=0.1)
    grid.grid_init()
    # set up x,y,z of grid
    x = grid.x
    y = grid.y
    z = grid.z
    
    # calculate densities
    mo_list, mo_info = extras.calc_mo(
    read_file,orbital_string_input,ofid='vis/mo')
    
    return x, y, z, mo_list, mo_info
    
    
    

