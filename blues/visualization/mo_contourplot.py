from . import *
from choose_calculate_mo import *
from visualize_mo_x3d import *
import matplotlib.pyplot as plt


def make_mesh_grid(selected_mo, filename=None, before_homo=None, 
after_homo=None, extend=7.0, output_name='charges', comments=''):
"""
Parameters
----------

Returns
----------

"""
    
    if mo_list 
    x, y, z, mo_list, mo_info = calculate_densities(filename, 
    before_homo, after_homo, extend=7.0, output_name='charges', comments='')
    
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
    
    return XZ, YZ, XY, YY, XX, YX, mo_list_xy, mo_list_xz, mo_list_yz
    
def contour_plots(contour_lines=10, linewidth=0.5, cmap='seismic'):
"""
Parameters
----------

Returns
----------

"""
    XZ, YZ, XY, YY, XX, YX, mo_list_xy, mo_list_xz, mo_list_yz
    = make_mesh_grid()
    
    # set up grid of 3x1 contour plots of electron density for selected MO
    f, (pic1, pic2, pic3) = 
    plt.subplots(3,1,sharex=True,sharey=True,figsize=(6,16))
    
    # figure 1
    pic1.contour(XZ,YZ,mo_list_xy,contour_lines,linewidths=0.5,colors='k')
    pic1.contourf(XZ,YZ,mo_list_xy,contour_lines,
    cmap=cmap,vmax=abs(mo_list_xy).max(),vmin=-abs(mo_list_xy).max())
    pic1.set_xlabel('x',fontsize=16)
    pic1.set_ylabel('y',fontsize=16)
    
    # set title for whole figure
    pic1.set_title('MO is {0}. MO energy is {1}'.
                format((orbital_list[selected_mo]),(mo_info['mo_spec'][selected_mo]['energy']),fontsize=18))
    
    # figure 2
    pic2.contour(XY,YY,mo_list_xz,contour_lines,linewidths=0.5,colors='k')
    pic2.contourf(XY,YY,mo_list_xz,contour_lines,
    cmap='seismic',vmax=abs(mo_list_xz).max(),vmin=-abs(mo_list_xz).max()) 
    pic2.set_xlabel('x',fontsize=16)
    pic2.set_ylabel('z',fontsize=16)
    
    # figure 3
    pic3.contour(XX,YX,mo_list_yz,contour_lines,linewidths=0.5,colors='k')
    pic3.contourf(XX,YX,mo_list_yz,contour_lines,
    cmap='seismic',vmax=abs(mo_list_yz).max(),vmin=-abs(mo_list_yz).max())
    pic3.set_xlabel('y',fontsize=16)
    pic3.set_ylabel('z',fontsize=16)     
    
    # plot
    return f.show()
    
    
    