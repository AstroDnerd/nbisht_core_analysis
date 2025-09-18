
from starter2 import *
from collections import defaultdict
import scipy
import colors
import track_loader as TL
import hair_dryer
reload(hair_dryer)


def more_fractal(this_looper,core_list=None, external_ax=None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    if external_ax is None:
        fig,ax=plt.subplots(1,1)
    else:
        ax = external_ax
    x_ext=extents()
    y_ext=extents()
    z_ext=extents()
    mini_scrubbers={}
    times = thtr.times
    times.shape = times.size,1
    dont_axbonk=False
    for core_id in core_list:
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
        ms.particle_pos(core_id)

        mini_scrubbers[core_id]=ms
        x_ext( ms.particle_x)
        y_ext( ms.particle_y)
        z_ext( ms.particle_z)
    ds=this_looper.load(0)
    for LOS in [2]:
        x = [1,0,1][LOS] # Using [1,0,1] and [2,2,0] 
        y = [2,2,0][LOS] # unfolds nicely.

        for ncore,core_id in enumerate(core_list):
            print('a')
            dx_0 = 1/128
            ms=mini_scrubbers[core_id]
            left  = nar([x_ext.minmax[0],y_ext.minmax[0],z_ext.minmax[0]])
            right = nar([x_ext.minmax[1],y_ext.minmax[1],z_ext.minmax[1]])
            #quantize
            left = np.floor(left/dx_0)*dx_0
            right= np.ceil(right/dx_0+1)*dx_0
            nzones = np.floor((right-left)/dx_0).astype('int')

            p = [ms.particle_x.transpose(),ms.particle_y.transpose(),ms.particle_z.transpose()]
            cg = ds.covering_grid(0,left,nzones)
            sl=slice(None)

            fig,axes=plt.subplots(2,2,figsize=(12,12))
            ax0=axes[0][0];ax1=axes[0][1];ax2=axes[1][0];ax3=axes[1][1]

            thax=2
            theH=0;theV=1
            theXx='xyz'[theH]; theYx='xyz'[theV]; theZx='xyz'[thax]

            z_coord_list = np.sort(np.unique(cg[theZx])).v
            qfull = cg['density']
            XX1= np.unique(cg[theXx].mean(axis=thax))
            YY1= np.unique(cg[theYx].mean(axis=thax))
            X2,Y2=np.meshgrid(XX1,YY1)

            figA,axA=plt.subplots(1,1,figsize=(12,12))
            for nz,the_z in enumerate(z_coord_list):
                print('slicer',nz)
                axA.clear()
                z_coord = np.argmin( np.abs(z_coord_list-the_z))
                sl = [slice(None),slice(None),slice(None)]
                sl[thax]=z_coord
                sl=tuple(sl)
                qslice=qfull[sl].transpose()
                axA.pcolormesh(X2,Y2,qslice,shading='nearest', cmap='Reds')

                ok = (np.abs(p[thax][0,:] - the_z)<1/128)
                x1=p[theH][0,:][ok]
                y1=p[theV][0,:][ok]
                axA.scatter(x1,y1)
                figA.savefig('plots_to_sort/zslices_%s_c%04d_%04d'%(this_looper.sim_name,core_id,nz))


sims=['u501','u502','u503']

sims=['u502']
TL.load_tracks(sims)
import find_other_cores

for sim in sims:
        corelist = [214,74]
        more_fractal(TL.loops[sim], core_list=corelist)

