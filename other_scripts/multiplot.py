import matplotlib.pyplot as plt
import matplotlib
from yt.visualization._mpl_imports import FigureCanvasAgg
import numpy as na
import pdb
def mp(nx,ny,colorbar_size = 2, size_x=None,size_y=None, ylabel_size=1.5, right_margin=0,xlabel_size=1,title_size=0.5,dpi=300,
      debug = 0):
    """Make *nx* by *ny* plots of size *size_x* by *size_y*, that share common title, xlabel, ylabel, and colorbar. 
    ylabel is for the y axis, whose real estate comes out of the x dimension.  Vice-versa for x."""


    total_x = 1.0*(nx*size_x + ylabel_size + colorbar_size + right_margin)
    total_y = 1.0*(ny*size_y + xlabel_size + title_size)
    if size_x is None:
        size_x = 4
    if size_y is None:
        size_y = 4
    fr_x = size_x / total_x
    fr_y = size_y / total_y
    if debug > 0:
        print "Total", total_x, total_y, "pixels", total_x*dpi,total_y*dpi

    fig = matplotlib.figure.Figure((total_x,total_y), dpi=dpi)
    fig.set_canvas(FigureCanvasAgg(fig))
    axes_list = []
    counter=0
    for j in range(ny-1,-1,-1):
        counter += 1
        axes_list.append([])
        for i in range(nx):
            left = (ylabel_size+i*size_x)/total_x
            bottom = (xlabel_size + j * size_y)/ total_y
            ax_extent = [left,bottom,size_x/total_x,size_y/total_y]
            ax = fig.add_axes(ax_extent)
            axes_list[-1].append(ax)
            if debug >0:
                print "new ax: (l,b,r,t) (%f, %f, %f, %f)"%tuple(ax_extent)
                print "new axb: (x0,x1) (%f %f)"%(ax.bbox.x0,ax.bbox.x1)
    cbar = None
    if colorbar_size > 0:
        left = (ylabel_size + nx*size_x)/total_x
        bottom = 1.0*xlabel_size/total_y
        cbar_extents = [left,bottom,(colorbar_size-ylabel_size)/total_x, ny*size_y/total_y]
        if debug > 0:
            print cbar_extents
        cbar = fig.add_axes( cbar_extents)
    return fig,axes_list,cbar
            


