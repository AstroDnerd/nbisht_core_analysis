from starter2 import *
import matplotlib.image as mpimg

import data_locations as dl
reload(dl)
plt.close('all')

file_list=glob.glob('%s/*h5'%dl.sixteen_frame)
out_prefix = 'CHANGE_THIS_PREFIX'

out_prefix='tracks'


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
xext = extents()
yext = extents()
axis = 'x'
do_art=False
if 1:
    used_cores=[]
    fig,axes=plt.subplots(1,1,figsize=(8,8))
    ax=axes
    ax.plot([0,1,1,0,0],[0,0,1,1,0],c=[0.5]*3)


    for core_id in core_list:
        plt.close('all')
        ms = trackage.mini_scrubber(thtr,core_id)
        if ms.nparticles < 3:
            continue
        if axis == 'x':
            the_x = ms.mean_y
            the_y = ms.mean_z
            xlabel='y'
            ylabel='z'
        if axis == 'y':
            the_x = ms.mean_z
            the_y = ms.mean_x
            xlabel='z'
            ylabel='x'
        if axis == 'z':
            the_x = ms.mean_x
            the_y = ms.mean_y
            xlabel='x'
            ylabel='y'
        ax.plot( the_x, the_y,c='k',linewidth=0.1)
        ax.scatter( the_x[0], the_y[0],c='k',s=1)
        xext(ms.mean_y)
        yext(ms.mean_z) 



    delta=0.1   
    delta = max([np.abs(xext.minmax[0]), np.abs(1-xext.minmax[1]),np.abs(yext.minmax[0]), np.abs(1-yext.minmax[1])])
    axbonk(ax,xlabel=xlabel,ylabel=ylabel,xlim=[-delta,1+delta],ylim=[-delta,1+delta])
    if do_art:
        axbonk(ax,xlabel='',ylabel='')
        ax.set_xticks([])
        ax.set_yticks([])
    ax.set_aspect('equal')
    fig.savefig('plots_to_sort/centroids_proj_%s.pdf'%axis)


    
