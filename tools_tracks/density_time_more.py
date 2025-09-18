from starter2 import *
import matplotlib.image as mpimg

import data_locations as dl
reload(dl)
plt.close('all')

file_list=glob.glob('%s/*h5'%dl.sixteen_frame)
out_prefix = 'CHANGE_THIS_PREFIX'

out_prefix='test'
#file_list = file_list[:2]


if 'this_looper' not in dir():
    this_looper=looper.core_looper(directory=dl.enzo_directory)
    for nfile,fname in enumerate(file_list):
        this_looper.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    thtr = this_looper.tr
    thtr.sort_time()
all_cores = np.unique(thtr.core_ids)
rm = rainbow_map(len(all_cores))
tsorted = thtr.times[asort]

core_list = all_cores
fast_cores = []
#core_list = [0,1,25,79]
FallThreshold = 5e3
for core_id in core_list:
    n0=0

    fig,ax=plt.subplots(1,1)
    ms = trackage.mini_scrubber(thtr,core_id)
    density = thtr.c([core_id],'density')

    tmap=rainbow_map(ms.ntimes)
    norm = mpl.colors.Normalize()
    norm.autoscale( np.log10(density[:,n0]))
    cmap = mpl.cm.jet
    color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)

    first_time_index=-1
    if density.max() > FallThreshold:
        whr = np.where(density > FallThreshold )
        first_time_index = whr[1].min()

    for npart in list(range(ms.nparticles)):
        c = color_map.to_rgba(np.log10(density[npart,n0]))
        ax.plot( thtr.times, density[npart,:],c=c,linewidth=.1)#linestyle=':')

    if first_time_index>0:  
        ax.scatter(thtr.times[first_time_index], density[:,first_time_index].max(),marker='*')

    mean_densities = density.mean(axis=0)[:]
    ax.plot(tsorted, mean_densities,c='k')

    t0 = thtr.times[0]

    t1 = thtr.times[-1]
    rho0 =1.1 
    rho0 = mean_densities[0]


    rho1 = density.max() 
    alpha = 1.8
    tc =t1*(1-(rho1/rho0)**(-1./alpha))**-0.5
    #G=np.pi*4#1620./(4*np.pi)
    G = 1620/(4*np.pi)
    tff_global = np.sqrt(3*np.pi/(32*G*1))
    tff_local = np.sqrt(3*np.pi/(32*G*rho0))
    timescale = tff_local
    ax.plot([tff_local]*2,[0.1,1e6],c='r')
    rhot = rho0*(1-(tsorted/timescale)**2)**-alpha
    rho_c = 3*np.pi/(32*G*tc**2)

    if thtr.times[first_time_index] < tff_local:
        ax.set_title("THIS ONE")
        fast_cores.append(core_id)
    else:
        ax.set_title("no")

    ok = np.isnan(rhot)==False
    ax.plot( tsorted[ok], rhot[ok], c='r',label=r'$tc/tff = %0.2e$'%(tc/tff_local))
    ax.legend(loc=0)
    axbonk(ax,xscale='linear',yscale='log',xlabel='t',ylabel=r'$\rho$')
    oname = "%s/%s_density_5_c%04d"%(dl.output_directory,out_prefix,core_id)
    fig.savefig(oname)
    print("Saved "+oname)
    plt.close(fig)
