from starter2 import *
import data_locations as dl
import davetools
reload(looper)
reload(trackage)
reload(davetools)

plt.close('all')

ratarray=[]
file_list=glob.glob('%s/*h5'%dl.sixteen_frame)[:10]
file_list = glob.glob('../Datasets/all_primitives/all_primitives_c000[01].h5')
file_list = glob.glob('../Datasets/all_primitives/all_primitives_c*.h5')
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

plt.close('all')
fig, axd1=plt.subplots(1,1)
fig4, ax4=plt.subplots(2,2)
ax40 = ax4[0][0]
ax41 = ax4[0][1]
ax42 = ax4[1][0]
ax43 = ax4[1][1]

do_vel_extent=False
if 'vel_ext' not in dir():
    vel_ext = extents()
    do_vel_extent=True

for nc,core_id in enumerate(core_list):

    #miniscrubber computes distance, r^2, several other quantities
    ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)

    if do_vel_extent:
        vel_ext(ms.rel_vx)
        vel_ext(ms.rel_vy)
        vel_ext(ms.rel_vz)
        vel_ext(ms.rel_v2**0.5)
        continue

    for ax in [axd1,ax40,ax41,ax42,ax43]:
        ax.clear()    


    tmap=rainbow_map(ms.ntimes)
    if ms.nparticles == 1:
        continue

    asort =  np.argsort(thtr.times)
    density = thtr.c([core_id],'density')
    if (asort != sorted(asort)).any():
        print("Warning: times not sorted.")
    n0=asort[0]
    tsorted = thtr.times[asort]


    for n_count,n_time in enumerate(asort):
        time=thtr.times[n_time]
        if time == 0:
            continue
        c=tmap(n_count,ms.nparticles)
        this_r=ms.r[:,n_time]+0
        this_r[ this_r < 1./2048] = 1./2048
        r_un = nar(sorted(np.unique(this_r)))

        if 0:
            axd1.scatter(ms.vr_rel,ms.rel_v2**0.5,c=c[0:1],label=thtr.times[n_time],s=0.1)
        if  1:
            outname4 = '%s/vi_r_symlog_c%04d'%(dl.output_directory,core_id)
            ax40.scatter(this_r, ms.rel_vx[:,n_time], c=c,s=0.1)
            ax41.scatter(this_r, ms.rel_vy[:,n_time], c=c,s=0.1)
            ax42.scatter(this_r, ms.rel_vz[:,n_time], c=c,s=0.1)
            ax43.scatter(this_r, ms.rel_v2[:,n_time]**0.5, c=c,s=0.1)
        if  0:
            #!!!!
            outname4 = '%s/vi_r_raw_c%04d'%(dl.output_directory,core_id)
            ax40.scatter(this_r, ms.raw_vx[:,n_time], c=c,s=0.1)
            ax41.scatter(this_r, ms.raw_vy[:,n_time], c=c,s=0.1)
            ax42.scatter(this_r, ms.raw_vz[:,n_time], c=c,s=0.1)
            ax43.scatter(this_r, ms.raw_v2[:,n_time]**0.5, c=c,s=0.1)

        if 0:
            v_radial_average = np.mean( ms.vr_rel[:,n_time])
            rel_vx_mean = np.mean(ms.rel_vx[:,n_time])
            rel_vy_mean = np.mean(ms.rel_vy[:,n_time])
            rel_vz_mean = np.mean(ms.rel_vz[:,n_time])
            v_total_average = np.sqrt(rel_vx_mean**2+rel_vy_mean**2+rel_vz_mean**2)
            axd1.scatter(v_radial_average, v_total_average, marker='s',c='k')
            #axd1.plot([-15,15],[1,1],c='k')
            #axd1.plot([-15,15],[-1,-1],c='k')
            axd1.plot([-10,0],[10,0],c='k')
            axd1.plot([0,10],[0,10],c='k')
        if 0:
            #axd1.scatter(ms.rel_v2[:,n_time],ms.vr_rel[:,n_time],c=c,label=thtr.times[n_time],s=0.1)
            abs_vr = np.abs(ms.vr_rel[:,n_time])
            #axd1.scatter(this_r,abs_vr,c=c,label=thtr.times[n_time],s=0.1)
            axd1.scatter(ms.vr_rel, ms.vr_rel/ms.rel_v2**0.5,c=c[0:1],label=thtr.times[n_time],s=0.1)
            harmonic_r =  10**( np.mean( np.log10(this_r)))
            mean_vr = np.mean(abs_vr)
            

            v_radial_average = np.mean( ms.vr_rel[:,n_time])
            rel_vx_mean = np.mean(ms.rel_vx[:,n_time])
            rel_vy_mean = np.mean(ms.rel_vy[:,n_time])
            rel_vz_mean = np.mean(ms.rel_vz[:,n_time])
            v_total_average = np.sqrt(rel_vx_mean**2+rel_vy_mean**2+rel_vz_mean**2)
            print(v_total_average)
            rat = v_radial_average/v_total_average
            if rat<-1:
                pdb.set_trace()
            ratarray.append(rat)
            axd1.scatter(v_radial_average, rat, marker='s',c='k')
            #axd1.plot([-10,0],[10,0],c='k')
            #axd1.scatter( harmonic_r, mean_vr, c=c[0:1], marker='*')

    limits = np.abs(vel_ext.minmax).max()
    limits = [-limits,limits]
    #axd1.plot([ms.r.min(),ms.r.max()], [1,1], c=[0.5]*4)
    labs = ['vx','vy','vz','vtotal']
    for iii,ax4i in enumerate([ax40,ax41,ax42,ax43]):
        davetools.axbonk(ax4i,xlabel='r',ylabel=labs[iii],ylim=limits)
        ax4i.set_xscale('symlog',linthreshx=1./2048)
        ax4i.set_yscale('symlog',linthreshy=2)
        ax4i.set_xlim([0,1])
    ax43.set_ylim([0,limits[1]])
    if 0:
        davetools.axbonk(axd1,xscale='linear',yscale='linear',xlabel='vr_rel',ylabel=r'$v_r/v_{total}$')
                         #xlim=[-15,15], ylim=[0,15])
        outname = '%s/ratio_c%04d'%(dl.output_directory,core_id)
    if 0:
        davetools.axbonk(axd1,xscale='linear',yscale='linear',xlabel='vr_rel',ylabel=r'$rel_v2$',
                         xlim=[-15,15], ylim=[0,15])
        outname = '%s/vr_vtotal_c%04d'%(dl.output_directory,core_id)
    if 0:
        davetools.axbonk(axd1,xscale='log',yscale='log',xlabel='r',ylabel=r'$v$',
                         xlim=r_extents.minmax, ylim=rho_extents.minmax)
        axd1.set_xscale('symlog',linthreshx=1./2048)
        axd1.set_xlim([0,1])
        axd1.set_yscale('symlog',linthreshy=1.)
        axd1.set_ylim([0,15])
        outname = '%s/vr_vtotal_c%04d'%(dl.output_directory,core_id)

    #outname = '%s/vr_vtotal_c%04d'%(dl.output_directory,core_id)
    #fig.savefig(outname)
    #print("saved "+outname)
    plt.close(fig)

    fig4.savefig(outname4)
    print("saved "+outname4)
    plt.close(fig4)
