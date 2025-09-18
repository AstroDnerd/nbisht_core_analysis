from starter2 import *
import data_locations as dl
import davetools
reload(davetools)

plt.close('all')


file_list=glob.glob('%s/*h5'%dl.sixteen_frame)
file_list=glob.glob('../Datasets/all_primitives/*h5')
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
    field_extents=davetools.extents()
    r_extents=davetools.extents()
    for nc,core_id in enumerate(all_cores):
        ms = trackage.mini_scrubber(thtr,core_id)
        if ms.nparticles == 1:
            continue
        density = thtr.c([core_id],'density')
        field = thtr.c([core_id],'magnetic_field_strength')
        field_extents(field)
        rho_extents(density)
        r_extents(ms.r)
theta_extents=davetools.extents()
alphab_list=[]
for nc,core_id in enumerate(core_list):

    ms = trackage.mini_scrubber(thtr,core_id)
    if ms.nparticles == 1:
        continue

    field = thtr.c([core_id],'magnetic_field_strength')
    bx = thtr.c([core_id],'magnetic_field_x')
    by = thtr.c([core_id],'magnetic_field_y')
    bz = thtr.c([core_id],'magnetic_field_z')
    vx = thtr.c([core_id],'velocity_x')
    vy = thtr.c([core_id],'velocity_y')
    vz = thtr.c([core_id],'velocity_z')

    bb =np.sqrt(bx*bx+by*by+bz*bz)
    vv =np.sqrt(vx*vx+vy*vy+vz*vz)
    BdotV = bx*vx+by*vy+bz*vz
    costheta=BdotV/(bb*vv)
    theta_extents(costheta)

    #miniscrubber computes distance, r^2, several other quantities
    tmap=rainbow_map(ms.ntimes)

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

        axd1.scatter(density[:,n_time],np.abs(costheta[:,n_time]),c=c,label=thtr.times[n_time],s=0.1)
        #axd1.plot(r_un, 100*(r_un/1e-2)**-2,c='k',linewidth=0.1)

    davetools.axbonk(axd1,xscale='log',yscale='linear',ylabel=r'$cos \theta$',xlabel=r'$\rho$',
                     ylim=[0,1], xlim=rho_extents.minmax)
    outname = '%s/density_abs_angle_c%04d'%(dl.output_directory,core_id)
    fig.savefig(outname)
    print("saved "+outname)
    plt.close(fig)

