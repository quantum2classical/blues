import numpy as np
import pandas as pd
import os

try:
    import orbkit as ok
except ImportError:
    import sys
    sys.path.append('C:\orbkit-cython')
from orbkit import read, grid, atomic_populations, options, main, extras
from orbkit.output import pdb_creator ,xyz_creator
from orbkit.display import init_display, display

# set up global 
global density_done, atoms_done
density_done = False
atoms_done = False

# this submodule should output:
# 1. strings of the number of MOs
# 2. a list of which MOs will be calculated
# 3. calculated density plots for all MOs
# 4. 3d grid with x, y, and z values

def read_molden(filename, itype='molden', all_mo=True):
    """
    parameters
    ----------
    filename: string of molden file to be read
    
    itype: string of file type, should be molden
    
    all_mo: bool, read all MOs (it's a quick calculation, should be True)
    
    returns
    ----------
    read_file: qc class object of read molden file
    
    Also prints which MO HOMO is. Example: 
    'HOMO is MO 50 out of 100'
    """
    
    # make sure right file type and itype matches   
    assert isinstance(filename, str), 'filename must be molden file string'
    assert filename[-6:-1]+filename[-1] == itype,'file type must match itype argument'
    
    # read molden file
    read_file = read.main_read(filename, itype='molden', all_mo=True)
    #find out how many MOs and where HOMO is by finding when MO occupancy = 0
    for x in range(len(read_file.mo_spec)):
        if read_file.mo_spec[x]['occ_num'] != 0.0:
            continue
        else:
            where_is_homo = 'HOMO is MO number ' + str(x - 1) + \
              ' out of ' + str(len(read_file.mo_spec))
            break
            
    print(where_is_homo)
    
    return read_file


def my_mo_list(before_homo, after_homo):
    """
    parameters
    ----------
    filename: string of molden file to be read
    
    before_homo: integer, how many MOs to calculate before HOMO. Can be < 0 
    if only MOs after homo are to be calculated. Must stay within MO range
    
    after_homo: integer, how many MOs to caluclate after HOMO

    
    returns:
    ----------
    
    orbital_string_input: string input used in calc_mo function
    
    Also prints list of strings of MOs which the user has chosen to be 
    calculated
    """
    
    # ensure correct input type
    assert isinstance(before_homo, int), 'input must be integer'
    assert isinstance(after_homo, int), 'input must be integer'
    

    if before_homo and after_homo < 0:
        raise Exception('list must begin at lowest MO, this list would be reversed')    
    elif after_homo < 0 and abs(after_homo) > abs(before_homo):
        raise Exception('if after_homo is negative, absolute value must be less or equal to before_homo')
    elif before_homo < 0 and abs(before_homo) > abs(after_homo):
        raise Exception('if before_homo is negative, absolute value must be less or equal to after_homo')
    else:
        pass
        
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
    
    return orbital_list, orbital_string_input



def get_atom_positions(filename, comments=''):
    """
    parameters
    ----------
    filename: string of .molden file to be read
    
    comments: string of comments embedded in pdb file
        
    returns
    ----------
    xyz_c, xyz_o, xyz_n: np.array of parsed atomic positions from pdb 
    file. One array for N, O, and the rest C and any other atom
    
    """
    
    # output name will be filename_charges.pdb
    output_name = filename[0:len(filename)-7] + str('_atom_positions')
    # read molden file
    read_file = read_molden(filename, itype='molden', all_mo=True) 
    
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
    atoms=atoms.loc[3:len(atoms)-3]
    atoms=atoms.iloc[:,2:6]
    # remove hydrogens
    atoms = [atoms.atom != 'H']
    
    # unit coversion
    ang2au = 0.529
    # scale to make units a.u.
    atoms['x']=atoms['x']/ang2au
    atoms['y']=atoms['y']/ang2au
    atoms['z']=atoms['z']/ang2au
    
    # create individual dataframes for O, N atoms
    atoms_o = atoms[atoms.atom=='O']
    atoms_n = atoms[atoms.atom =='N']
    # remove O,N from original dataframe
    atoms = atoms[atoms.atom != 'O']
    atoms = atoms[atoms.atom !='N'] 
    
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
    
    # move file to atom_positions directory
    os.rename(output_name+'.pdb', 'atom_positions/'+output_name+'.pdb')
    
    atoms_done == True
    
    return xyz_c, xyz_o, xyz_n
    

def calculate_densities(filename, before_homo, after_homo, extend=2.0,
                        numproc=None):
    """
    parameters
    ---------
    filename: string of molden file to be read
    
    before_homo: integer, how many MOs to calculate before HOMO. Can be < 0 
    if only MOs after homo are to be calculated. Must stay within MO range
    
    after_homo: integer, how many MOs to caluclate after HOMO
    
    extend: float, how much to extend grid so MOs can be visualized. 
    2.0 default. Warning - extending by large amount will drastically slow down
    calculation
    
    numproc: int, number of processors available. Default is None
    
    returns
    ---------
    x, y, z: np.arrays making up 3d grid on which densities will be plotted
    
    mo_list: list of dictionaries containing MO information, including densities
    
    mo_info: dictionary of info for each MO
    
    also prints: 
    where_is_homo -  number of MOs and where HOMO lies
    orbital_list: list of MOs that user chose to calculate
   
    """
    
    orbital_list, orbital_string_input = my_mo_list(filename, before_homo, after_homo)
    read_file = read_molden(filename, itype='molden', all_mo=True)
    
    # create file with atom positions
    
    # adjust grid to fit molecule + extended units
    # this prevents calc_mo function from automatically setting up grid     
    grid.adjust_to_geo(read_file,extend=extend,step=0.1)
    grid.grid_init()
    # set up x,y,z of grid
    x = grid.x
    y = grid.y
    z = grid.z
    
    # calculate densities
    mo_list, mo_info = extras.calc_mo(
    read_file, orbital_string_input, numproc=numproc)
    
    density_done == True
    return x, y, z, mo_list, mo_info
    
    
    

