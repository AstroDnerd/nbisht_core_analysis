from starter1 import *
import davetools
def simple_phase(field1,field2,log=False,ax=None,nBins=64):
    ext=davetools.extents()
    ext(field1)
    ext(field2)
    if log:
        bins=np.geomspace(ext.minmax[0],ext.minmax[1],nBins)
    else:
        bins=np.linspace(ext.minmax[0],ext.minmax[1],nBins)
    hist,xbins,ybins=np.histogram2d( field1,field2, bins=[bins,bins])
    helper( hist, xbins, ybins,ax=ax)
def helper(h_in,xbins_in,ybins_in, cmap_name = 'viridis', zlim=None, ax=None,transpose=False, **pcolormesh_args):
    #takes the output of np.histogram2d (or any other 2d histogram)
    #xbins is 1 larger than h.size[0].

    xbins = xbins_in+0
    ybins = ybins_in+0
    h=h_in+0

    if transpose:
        h = h.transpose()
        temp = xbins.transpose()
        xbins=ybins.transpose()
        ybins=temp


    if xbins.size > h.shape[0]:
        xbins = 0.5*(xbins[1:] + xbins[:-1])
        ybins = 0.5*(ybins[1:] + ybins[:-1])

    nx = len(xbins) ; ny=len(ybins)
    TheX = np.r_[(ny)*[xbins]].transpose()
    TheY = np.r_[(nx)*[ybins]]

    if zlim is None:
        zmin = h[h>0].min()
        zmax = h.max()
    else:
        zmin = zlim[0]
        zmax = zlim[1]
    norm = mpl.colors.LogNorm( vmin =zmin, vmax=zmax)
    cmap = copy.copy(mpl.cm.get_cmap(cmap_name))
    cmap.set_under('w')
    if 'shading' not in pcolormesh_args:
        pcolormesh_args['shading'] = 'nearest'
    pcolormesh_args['norm']=norm
    pcolormesh_args['cmap']=cmap
    if ax is not None:

        ploot=ax.pcolormesh( TheX, TheY, h, **pcolormesh_args)

    output = {'TheX':TheX, 'TheY':TheY, 'norm':norm,'cmap':cmap, 'plot':ploot}
    return output




