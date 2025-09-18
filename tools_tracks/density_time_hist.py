import matplotlib.colors as colors

from starter2 import *
import data_locations as dl
reload(dl)
plt.close('all')

file_list=glob.glob('%s/*h5'%dl.sixteen_frame)
out_prefix = 'CHANGE_THIS_PREFIX'

out_prefix='test'
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

core_list = all_cores
fig, ax= plt.subplots(1,1)
nbins=100
for core_id in core_list:

    ax.clear()
    thhist = np.zeros([nbins,len(thtr.times)])
    ms = trackage.mini_scrubber(thtr,core_id)
    density = thtr.c([core_id],'density')
    tmap=rainbow_map(ms.ntimes)
    asort =  np.argsort(thtr.times)
    mind = density.min()
    maxd = density.max()
    bins = np.logspace(np.log10(mind), np.log10(maxd),nbins+1)
    fig,ax = plt.subplots(1,1)
    for itime,ntime in enumerate(asort):
        thhist[:,itime]=np.histogram(density[:,ntime],bins)[0]
    norm = colors.LogNorm(vmin=1,vmax=thhist.max())
    ploot=ax.imshow(thhist,origin='lower',interpolation='nearest',norm=norm)
    ax.set_aspect(.3)
    cbar = plt.colorbar(ploot,norm=norm)
    cbar.set_clim(vmin=1,vmax=thhist.max())
    cbar.cmap.set_under('w')

    outname = '%s/density_time_hist_c%04d.png'%(dl.output_directory,core_id)
    fig.savefig(outname )
    print("wrote "+outname)
if 0:
    fig2,ax2=plt.subplots(1,1)
    nparticles = density.shape[0]
    times_ordered=thtr.times[asort]
    for nt,t in enumerate(times_ordered):
        ax2.scatter([t]*nparticles,density[:,asort[nt]],c='k',s=.01)
    ax2.set_yscale('log')
    ax2.set_ylim(0.01,1e8)
    fig2.savefig('density_time_scatter.png')
        
plt.close('all')
