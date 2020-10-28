"""

Computing moments of inertia for p19 cores.

This is not a wise thing to do, as the cores are not
properly following the flow.  A full treatment of 
Virial Theorem including boudnary terms is necessray here



"""
from starter2 import *
import matplotlib.image as mpimg

import data_locations as dl
reload(dl)
plt.close('all')

file_list=glob.glob('%s/*h5'%dl.sixteen_frame)
out_prefix = 'CHANGE_THIS_PREFIX'

out_prefix='tracks'
file_list = file_list[:2]


if 'this_looper' not in dir():
    this_looper=looper.core_looper(directory=dl.enzo_directory)
    for nfile,fname in enumerate(file_list):
        this_looper.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    thtr = this_looper.tr
    thtr.sort_time()
all_cores = np.unique(thtr.core_ids)
rm = rainbow_map(len(all_cores))

tsorted = thtr.times
core_list = all_cores
moi_time = {}
mass_time = {}
if 1:
    used_cores=[]
    for core_id in core_list:
        plt.close('all')
        ms = trackage.mini_scrubber(thtr,core_id)
        if ms.nparticles < 3:
            continue
        density = thtr.c([core_id],'density')
        cell_volume = thtr.c([core_id],'cell_volume')
        used_cores.append(core_id)
        moi_time[core_id] = []
        mass_time[core_id] = []
        n0=0

        for nf, nframe in enumerate(thtr.frames):
            dx=1./2048
            nx=2048
            x =np.floor(thtr.c([core_id],'x')/dx)[:,nf]
            y =np.floor(thtr.c([core_id],'y')/dx)[:,nf]
            z =np.floor(thtr.c([core_id],'z')/dx)[:,nf]
            this_density = density[:,nf]
            cell_volume = thtr.c([core_id],'cell_volume')[:,nf]
            index = x + nx*(y * nx*z)
            ar = np.argsort(index)
            rs = np.argsort(ar)
            isorted=index[ar]
            mask = np.ones_like(this_density,dtype='bool')
            mask[1:] = isorted[1:]-isorted[:-1] != 0
            mask2 = mask[ rs]
            this_moi = (this_density[mask2]*cell_volume[mask2]*ms.r[:,nf][mask2]**2).sum()
            moi_time[core_id].append(this_moi)
            mass_time[core_id].append((this_density[mask2]*cell_volume[mask2]).sum())
        fig,axes=plt.subplots(1,2)
        ax=axes[0]
        axb=axes[1]

        ax.plot(thtr.times, moi_time[core_id] ,c='k')
        axb.plot(thtr.times, mass_time[core_id] ,c='k')
        axbonk(ax,xlabel=r'$t$',ylabel=r'$MOI$')
        axbonk(axb,xlabel=r'$t$',ylabel=r'$MOIb$')
        outname = 'plots_to_sort/moi_time_c%04d'%core_id
        fig.savefig(outname)
        print('saved '+outname)

