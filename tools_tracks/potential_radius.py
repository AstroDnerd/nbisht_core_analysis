from starter2 import *
import data_locations as dl
import davetools
reload(davetools)

plt.close('all')


file_list=glob.glob('%s/*h5'%dl.sixteen_frame)
#for debug purposes you may want a reduced list 
#file_list=file_list[:3]    

if 'this_looper' not in dir():
    this_looper=looper.core_looper(directory=dl.enzo_directory)
    for nfile,fname in enumerate(file_list):
        this_looper.load_loop(fname)
        print( "Reading file %d of %d"%(nfile,len(file_list)))
    thtr = this_looper.tr
    thtr.sort_time()
    all_cores = np.unique(thtr.core_ids)

core_list=all_cores
rm = rainbow_map(len(all_cores))

if 'rho_extents' not in dir():
    rho_extents=davetools.extents()
    r_extents=davetools.extents()
    for nc,core_id in enumerate(all_cores):
        ms = trackage.mini_scrubber(thtr,core_id)
        if ms.nparticles == 1:
            continue
        density = thtr.c([core_id],'density')
        rho_extents(density)
        r_extents(ms.r)

for nc,core_id in enumerate(core_list):

    #miniscrubber computes distance, r^2, several other quantities
    ms = trackage.mini_scrubber(thtr,core_id)
    tmap=rainbow_map(ms.ntimes)
    if ms.nparticles == 1:
        continue

    asort =  np.argsort(thtr.times)
    density = thtr.c([core_id],'density')
    if (asort != sorted(asort)).any():
        print("Warning: times not sorted.")
    n0=asort[0]
    tsorted = thtr.times[asort]

    fig, axd1=plt.subplots(1,1)

    for n_count,n_time in enumerate(asort):
        time=thtr.times[n_time]
        if time == 0:
            continue
        c=tmap(n_count,ms.nparticles)
        this_r=ms.r[:,n_time]+0
        r_un = nar(sorted(np.unique(this_r)))

        axd1.scatter(this_r,density[:,n_time],c=c,label=thtr.times[n_time],s=0.1)
        axd1.plot(r_un, 100*(r_un/1e-2)**-2,c='k',linewidth=0.1)

    davetools.axbonk(axd1,xscale='log',yscale='log',xlabel='r',ylabel=r'$\rho$',
                     xlim=r_extents.minmax, ylim=rho_extents.minmax)
    outname = '%s/density_radius_c%04d'%(dl.output_directory,core_id)
    fig.savefig(outname)
    print("saved "+outname)
    plt.close(fig)
