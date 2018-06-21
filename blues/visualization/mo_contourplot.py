from . import choose_calculate_mo
from choose_calculate_mo import *



def make_mesh_grid(selected_mo, filename=None, before_homo=None, 
after_homo=None, extend=7.0, output_name='charges', comments=''):
"""
Parameters
----------

Returns
----------

"""
    
    # check if density calculations have been done
    if density_done=False:
        x, y, z, mo_list, mo_info = calculate_densities(filename, 
        before_homo, after_homo, extend=7.0)
    else:
        pass
        
    # check if atomic positions have been found
    if mesh_done=False:
        get_atom_positions(filename, output_name='charges', comments='')
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
    
    return XZ, YZ, XY, YY, XX, YX, mo_slices
    
def contour_plots(contour_lines=10, linewidth=0.5, cmap='seismic'):
"""
Parameters
----------

Returns
----------

"""
    XZ, YZ, XY, YY, XX, YX, mo_slices = make_mesh_grid(....)
    
    # set up grid of 3x1 contour plots of electron density for selected MO
    f, (pic1, pic2, pic3) = 
    plt.subplots(3,1,sharex=True,sharey=True,figsize=(6,16))
    
    # figure 1
    pic1.contour(XZ,YZ,mo_slices[0],contour_lines,linewidths=0.5,colors='k')
    pic1.contourf(XZ,YZ,mo_list[0],contour_lines,
    cmap=cmap,vmax=abs(mo_slices[0]).max(),vmin=-abs(mo_slices[0]).max())
    pic1.set_xlabel('x',fontsize=16)
    pic1.set_ylabel('y',fontsize=16)
    
    # set title for whole figure
    pic1.set_title('MO is {0}. MO energy is {1}'.
                format((orbital_list[selected_mo]),(mo_info['mo_spec'][selected_mo]['energy']),fontsize=18))
    
    # figure 2
    pic2.contour(XY,YY,mo_slices[1],contour_lines,linewidths=0.5,colors='k')
    pic2.contourf(XY,YY,mo_slices[1],contour_lines,
    cmap='seismic',vmax=abs(mo_slices[1]).max(),vmin=-abs(mo_slices[1]).max()) 
    pic2.set_xlabel('x',fontsize=16)
    pic2.set_ylabel('z',fontsize=16)
    
    # figure 3
    pic3.contour(XX,YX,mo_slices[2],contour_lines,linewidths=0.5,colors='k')
    pic3.contourf(XX,YX,mo_slices[2],contour_lines,
    cmap='seismic',vmax=abs(mo_slices[2]).max(),vmin=-abs(mo_slices[2]).max())
    pic3.set_xlabel('y',fontsize=16)
    pic3.set_ylabel('z',fontsize=16)     
    
    # plot
    return f.show()
    
    
    