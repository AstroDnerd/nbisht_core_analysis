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
particle_slice=slice(None,None,10)
#core_list=all_cores[:5]
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
    tmap=rainbow_map(ms.ntimes)
    norm = mpl.colors.Normalize()
    norm.autoscale( np.log10(density[:,n0]))
    cmap = mpl.cm.jet
    color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)
    
    norm_theta = mpl.colors.Normalize()
    norm_theta.autoscale( [-1,1])
    cmap = mpl.cm.jet
    color_map_theta = mpl.cm.ScalarMappable(norm=norm_theta,cmap=cmap)

    fig, axd1=plt.subplots(1,1)
    ax=axd1

    for npart in list(range(ms.nparticles)[particle_slice]):
        #c = color_map.to_rgba(np.log10(density[npart,n0]))
        c1 = color_map_theta.to_rgba(costheta[npart,n0+1])
        the_y = costheta[npart,:][1:]
        the_t = thtr.times[1:]
        ax.scatter( the_t, the_y,c=[c1]*len(the_y),linewidth=.1)#linestyle=':')
        ax.plot( the_t, the_y,c=c,linewidth=.1)#linestyle=':')
    davetools.axbonk(axd1,xscale='linear',yscale='linear',ylabel=r'$cos \theta$',xlabel=r'$t$')
    outname = '%s/angle_time_c%04d'%(dl.output_directory,core_id)
    fig.savefig(outname)
    print("saved "+outname)
    plt.close(fig)

