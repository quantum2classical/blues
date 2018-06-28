from . import choose_calculate_mo
from choose_calculate_mo import *


def mo_to_visualize(selected_mo_string, contour1, contour2, notebook=True, 
                    opacity=0.3, scale_factor=0.5, wireframe=False,
                    saved_image='my_mo.png',filename=None, before_homo=None, 
                    after_homo=None ):
    """
    x3d rendering can only be used in a Jupyter Notebook for now. If the user 
    only wants a .png, change notebook=False
    
    Parameters
    ---------
    selected_mo: string, MO relative to HOMO user wants to view. Example, 
        'HOMO-30','LUMO', 'HOMO+6'. Must be in orbital list.
        
    contour1: float < 1.0 for isosurface of one phase. Contour values will be
    based on max value of density. 
    Example - contour = contour1 * mo_list[selected_mo].max()
    
    contour2: float < 1.0 for isofurface of opposite phase. Contour values will 
    be based on max value of density. 
    Example - contour = contour1 * mo_list[selected_mo].max()
    
    notebook: boolean, 
        
    Returns
    ---------
    1 x3d or png rendering with 2 isosurfaces and atomic positions
    """

    if density_done = False:
        calculate_densities(filename, before_homo, after_homo, extend=2.0)
    else:
        pass


    if notebook = True:
        # clear figure and re-initialize notebook with x3d for interactive,
        # run this each time you want a new figure    
        mlab.clf()
        mlab.init_notebook('x3d',500,500)
    else:
        pass
      
    read_file, orbital_string_input = 
    my_mo_list(filename, before_homo, after_homo)
    
    # convert string to index    
    selected_mo_string = selected_mo_string
    selected_mo = orbital_list.index(selected_mo_string)
    
    # contour spacing based on max density value
    contour1 = [mo_list[selected_mo].max()*contour]
    contour2 = [-mo_list[selected_mo].min()*contour]
    print(contour1, contour2)
    
    opacity = 0.3
    # auto-set the image boundaries from grid previously set
    extent= [x.min(),x.max(),y.min(),y.max(),z.min(),z.max()]

    # build pipeline
    src = mlab.pipeline.scalar_field(mo_list[selected_mo])
    
    isosurf = mlab.pipeline.iso_surface(
        src, contours=contour1, opacity=opacity, 
        color=(0, 0, 0.8), extent=extent)
    # isosurface of opposite phase
    isosurf = mlab.pipeline.iso_surface(
        src, contours=contour2, opacity=opacity, 
        color=(0.8, 0, 0), extent=extent)
    
    # plot atomic positions
    # make oxygens red
    isosurf = mlab.points3d(xo, yo, zo, color=(1.0,0.0,0.0), scale_factor=scale_factor)
    # make nitrogen blue
    isosurf = mlab.points3d(xn, yn, zn, color=(0.0,0.0,1.0), scale_factor=scale_factor)
    # all others (should be mostly carbon), white
    isosurf = mlab.points3d(xc, yc, zc, scale_factor=scale_factor)
        
    if wireframe=True:
        iso.actor.property.representation = 'wireframe'
    else:
        pass
    
        if notebook=False:
        mlab.axes()
        mlab.options.offscreen = True
        mlab.savefig(figure=isosurf, filename=selected_mo_string+'.png')
    else:
        pass
        
    return 