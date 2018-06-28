# -*- coding: utf-8 -*-
from . import choose_calculate_mo
from choose_calculate_mo import *



def make_mesh_grid(selected_mo, filename=None, before_homo=None, 
after_homo=None, extend=2.0, output_name='charges', comments=''):
    """
    Parameters
    ----------
    selected_mo: string, MO relative to HOMO user wants to view. Example, 
    'HOMO-30','LUMO', 'HOMO+6'. Must be in orbital list.
    
    filename
    Returns
    ----------
    XZ, YZ, XY, YY, XX, YX: numpy paired mesh grids for plotting the grid for
    contour lines
    
    mo_slices: list of 3 arrays, all 2D np.array slices of density data for
    plotting contour lines
    
    xyz: list containing 2D np.arrays of element positions for N, O, then C
    and all other elements    
    """
    
    # check if density calculations have been done
    if density_done = False:
        x, y, z, mo_list, mo_info = calculate_densities(filename, 
        before_homo, after_homo, extend)
    else:
        pass
        
    # check if atomic positions have been found
    if atoms_done=False:
        xyz_c, xyz_o, xyz_n =
        get_atom_positions(filename, output_name='charges', comments='')
        xyz = [xyz_c, xyz_o, xyz_n]
    else:
        pass
    
        
    # set up meshgrids to plot contours, make transpose so they match sizes
    XZ,YZ = np.meshgrid(x,y)
    XY,YY = np.meshgrid(x,z)
    XX,YX = np.meshgrid(y,z)
    XZ = XZ.T
    YZ = YZ.T
    XY = XY.T
    YY = YY.T
    XX = XX.T
    YX = YX.T
    # slice density data for contour plots
    mo_list_xy = mo_list[selected_mo][:,:,0]
    mo_list_xz = mo_list[selected_mo][:,0,:]
    mo_list_yz = mo_list[selected_mo][0,:,:]
    
    mo_slices = [mo_list_xy, mo_list_xz, mo_list_yz]
    
    return XZ, YZ, XY, YY, XX, YX, mo_slices, xyz
    
    
def contour_plots(selected_mo, contour_lines=10, linewidth=0.5, cmap='seismic', 
filename=None, before_homo=None, after_homo=None, extend=2.0, 
output_name='charges', comments=''):
    """
    Parameters
    ----------
    selected_mo: string, MO relative to HOMO user wants to view. Example, 
    'HOMO-30','LUMO', 'HOMO+6'. Must be in orbital list.
    
    contour_lines: integer, number of contour lines to divide min from max value.
    Default is 10
    
    linewidth: float, contour line widths. Default is 0.5
    
    cmap: string, any valid matplotlib colormap. Default is seismic.
    
    Returns
    ----------
    
    """
    # fontsize
    font = 16
    figsize=(6,16)
    
    # get mesh, density slices, and atom locations
    XZ, YZ, XY, YY, XX, YX, mo_slices, xyz = make_mesh_grid(selected_mo, 
    filename, before_homo, after_homo, extend, output_name='charges', comments='')
    # parse atom locations
    xyz[0] = xyz_c
    xyz[1] = xyz_o
    xyz[2] = xyz_n 
    
    # set up grid of 3x1 contour plots of electron density for selected MO
    f, (pic1, pic2, pic3) = 
    plt.subplots(3,1,sharex=True,sharey=True,figsize=figsize)
    
    # figure 1
    pic1.contour(XZ,YZ,mo_slices[0],contour_lines,linewidths=0.5,colors='k')
    pic1.contourf(XZ,YZ,mo_slices[0],contour_lines,
    cmap=cmap,vmax=abs(mo_slices[0]).max(),vmin=-abs(mo_slices[0]).max())
    pic1.set_xlabel('x',fontsize=16)
    pic1.set_ylabel('y',fontsize=16)
    
    # atom locations
    color=(0, 0, 0)
    pic1.scatter(xyz_c[0], xyz_c[1], marker=r"$ {} $".format('C'), c=color)
    pic1.scatter(xyz_n[0], xyz_n[1], marker=r"$ {} $".format('N'), c=color)
    pic1.scatter(xyz_o[0], xyz_o[1], marker=r"$ {} $".format('O'), c=color)
    
    pic1.set_xlabel('x(a.u.)',fontsize=font)
    pic1.set_ylabel('y(a.u.)',fontsize=font)
    
    # set title for whole figure
    pic1.set_title('Charge Density of {0}. MO energy is {1}'
    .format((orbital_list[selected_mo]),
    (mo_info['mo_spec'][selected_mo]['energy']),fontsize=18))
    
    # figure 2
    pic2.contour(XY,YY,mo_slices[1],contour_lines,linewidths=0.5,colors='k')
    pic2.contourf(XY,YY,mo_slices[1],contour_lines,
    cmap='seismic',vmax=abs(mo_list_xz).max(),vmin=-abs(mo_list_xz).max()) 
    
    # atom locations
    pic2.scatter(xyz_c[0], xyz_c[2], marker=r"$ {} $".format('C'), c=color)
    pic2.scatter(xyz_n[0], xyz_n[2], marker=r"$ {} $".format('N'), c=color)
    pic2.scatter(xyz_o[0], xyz_o[2], marker=r"$ {} $".format('O'), c=color)
    
    pic2.set_xlabel('x(a.u.)', fontsize=font)
    pic2.set_ylabel('z(a.u.)', fontsize=font)
    
    # figure 3
    pic3.contour(XX,YX,mo_slices[2],contour_lines,linewidths=0.5,colors='k')
    pic3.contourf(XX,YX,mo_slices[2],contour_lines,
    cmap='seismic',vmax=abs(mo_list_yz).max(),vmin=-abs(mo_list_yz).max())  
    # atom locations
    pic3.scatter(xyz_c[1], xyz_c[2], marker=r"$ {} $".format('C'), c=color)
    pic3.scatter(xyz_n[1], xyz_n[2], marker=r"$ {} $".format('C'), c=color)
    pic3.scatter(xyz_o[1], xyz_o[2], marker=r"$ {} $".format('C'), c=color)
    
    pic3.set_xlabel('y(a.u.)', fontsize=font)
    pic3.set_ylabel('z(a.u.)', fontsize=font)     
    
    # plot
    return f.show()
    
    