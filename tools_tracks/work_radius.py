from starter2 import *
import data_locations as dl
import davetools
reload(davetools)
from scipy import stats

plt.close('all')


#file_list=glob.glob('frame_40_c????.h5')
#file_list=glob.glob('frame_2_20_30*h5')
file_list=glob.glob('frame_40-125*h5')
file_list=glob.glob("../Datasets/SORTME/frame_40-125_c*.h5")[:30]
file_list=glob.glob("energy_test_full*h5")
#file_list=glob.glob('%s/track_indfix_sixteenframe_core_*.h5'%( "../Datasets/track_indfix_sixteenframe/"))
#for debug purposes you may want a reduced list 
#file_list=file_list[:3]    
output_prefix='Work1'

if 'this_looper' not in dir():
    this_looper=looper.core_looper(directory=dl.enzo_directory)
    for nfile,fname in enumerate(file_list):
        this_looper.load_loop(fname)
        print( "Reading file %s %d of %d"%(fname,nfile,len(file_list)))
    thtr = this_looper.tr
    thtr.sort_time()
    all_cores = np.unique(thtr.core_ids)

core_list=all_cores
rm = rainbow_map(len(all_cores))

if 1:
    if 'all_extents' not in dir():
        print("rerun grav extents")
        all_extents=davetools.extents()
        r_extents=davetools.extents()
        for nc,core_id in enumerate(all_cores):
            ms = trackage.mini_scrubber(thtr,core_id)
            if ms.nparticles == 1:
                continue
            grav_energy = thtr.c([core_id],'grav_energy')
            all_extents(grav_energy)
            r_extents(ms.r)

if 1:

    r_min = 0.5/2048 #   
    r_max = r_extents.minmax[1]*(1.00001)
    r_bins = np.linspace(0,r_max,64)
#r_bins = np.logspace(np.log(r_min),np.log(r_max),100)
    r_centers = 0.5*(r_bins[1:]+r_bins[:-1])
    #field_list = ['therm_energy','grav_energy','magnetic_energy','kinetic_energy']
    all_fields = ['momentum_flux','gravity_work','mag_work','pressure_work']
    field_list = ['momentum_flux','gravity_work','mag_work','pressure_work']
    #field_list = ['mag_work']
    labels={}
    density_extents=davetools.extents()
    density_extents(thtr.track_dict['density'])
    def get_field_color(field):
        color = 'rgbcmyk'[ all_fields.index(field)]
        return color

#field_list=['density']
#output_prefix='density'
    for nc,core_id in enumerate(core_list[:]):

        #miniscrubber computes distance, r^2, several other quantities
        ms = trackage.mini_scrubber(thtr,core_id)
        tmap=rainbow_map(ms.ntimes)
        if ms.nparticles == 1:
            continue

        asort =  np.argsort(thtr.times)

        if (asort != sorted(asort)).any():
            print("Warning: times not sorted.")
        n0=asort[0]
        tsorted = thtr.times[asort]

        fig, axd1=plt.subplots(1,1)

        frame_stuff = ""
        
        for n_count,n_time in enumerate(asort):
            time=thtr.times[n_time]
            if time == 0:
                continue
            cxx=tmap(n_count,ms.nparticles)
            this_r=ms.r[:,n_time]+0
            this_r[ this_r<1./2048] = 1./2048

            for nf,field in enumerate(field_list):
                c = get_field_color(field)
                #c = tmap(n_time)
                this_field = thtr.c([core_id],field)
                all_extents(this_field)
                arr = this_field[:,n_time]
                field_growth = arr >= 0
                field_death  = arr < 0
                if field_growth.sum() > 0:
                    bin_means_g, bin_edges_g, binnumber_g=stats.binned_statistic(this_r[field_growth], arr[field_growth], statistic='mean', bins=r_bins)
                    ok = ~np.isnan(bin_means_g)
                    print("OK %d %d"%(ok.sum(),ok.size))
                    label=None
                    if n_count == 0:
                        label=field
                    axd1.plot( r_centers[ok], bin_means_g[ok],c=c,label=label,linestyle='-')
                if field_death.sum() > 0:
                    bin_means_d, bin_edges_d, binnumber_d=stats.binned_statistic(this_r[field_death], arr[field_death], statistic='mean', bins=r_bins)
                    label=None
                    axd1.plot( r_centers, bin_means_d,c=c,label=label,linestyle='--')
                if True:
                    bin_means_b, bin_edges_b, binnumber_b=stats.binned_statistic(this_r, arr, statistic='mean', bins=r_bins)
                    print(bin_means_b)
                    label=None
                    if n_count == 0:
                        label=field
                    axd1.plot( r_centers, bin_means_b,c=c,label=label)

#                axd1.scatter(this_r,arr,c=cxx,label=labels.get(field,field),s=0.1,marker='*')
                print("N pos %d N neg %d Ntot %d FG %d"%( (arr>=0).sum(), (arr<0).sum(),arr.size, field_growth.sum()))

            #r_un = nar(sorted(np.unique(this_r)))
            #axd1.plot(r_un, 100*(r_un/1e-2)**-2,c='k',linewidth=0.1)

        #davetools.axbonk(axd1,xscale='log',yscale='log',xlabel='r',ylabel=r'$E_G,E_B$', xlim=r_extents.minmax, ylim=all_extents.minmax)
        davetools.axbonk(axd1,xscale='log',yscale='linear',xlabel='r',ylabel=r'$E_G,E_B$', xlim=[r_min,r_extents.minmax[1]], ylim=all_extents)
        #axd1.legend(loc=1)
        axd1.set_yscale('symlog',linthreshy=1e-3)
        axd1.set_ylim(all_extents.minmax)
        #axd1.legend(loc=0)
        outname = '%s/%s_radius_c%04d'%(dl.output_directory,output_prefix,core_id)
        axd1.set_title("core %d %s"%(core_id,frame_stuff))
        fig.savefig(outname)
        print("saved "+outname)
        plt.close(fig)
