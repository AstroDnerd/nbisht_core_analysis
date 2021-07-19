
clobber=False
import three_loopers as tl
import tools_tracks.density_time as DT
reload(DT)
import tools_tracks.neighbors as NE
reload(NE)
clobber=True
if 'density_tool1' not in dir():
    density_tool1=DT.trial(tl.looper1)
    density_tool1.run(do_all_plots=False)
    density_tool2=DT.trial(tl.looper2)
    density_tool2.run(do_all_plots=False)
    density_tool3=DT.trial(tl.looper3)
    density_tool3.run(do_all_plots=False)

if 'neighbor_tool1' not in dir() or clobber:
    neighbor_tool1=NE.neighbors(tl.looper1)
    neighbor_tool1.run()
if 'neighbor_tool2' not in dir() or clobber:
    neighbor_tool2=NE.neighbors(tl.looper2)
    neighbor_tool2.run()
if 'neighbor_tool3' not in dir() or clobber:
    neighbor_tool3=NE.neighbors(tl.looper3)
    neighbor_tool3.run()


if 1:
    fig,ax=plt.subplots(2,2)
    axex=ax.flatten()
    neighbors = [neighbor_tool1,neighbor_tool2,neighbor_tool3]
    densities = [density_tool1,density_tool2,density_tool3]
    for i in range(3):
        times = nar(densities[i].collapse_times)
        ok = times>0
        axex[i].scatter( neighbors[i].mean_density[ok], times[ok],label=neighbors[i].this_looper.out_prefix)
        axex[i].legend(loc=0)
        axbonk(axex[i],xlabel=r'$\rho_5$',ylabel=r'$t_{c}$',xscale='log',yscale='log')

    outname = "plots_to_sort/density_rate.pdf"
    fig.savefig(outname)
    print(outname)

if 1:
    def dump_core_vals(cores,vals,fname='out.h5',setname='values'):
        if len(cores) != len(vals):
            print("Core list and val list do not match")
            raise
        if os.path.exists(fname):
            mode = 'r+'
        else:
            mode = 'w'
        print("OPEN",fname, " mode ", mode)
        fptr=h5py.File(fname,mode)
        try:
            fptr.create_dataset("core_ids",data=cores)
            fptr.create_dataset(setname,data=vals)
        except:
            raise
        finally:
            fptr.close()
    for i in range(3):
        nt = neighbors[i]
        cores = nt.cores_used
        #dump_core_vals(cores, nt.mean_density, "%s_neighbors_mean_density.h5"%nt.this_looper.out_prefix,'mean_density')
        #dump_core_vals(cores, nt.NEXT_DISTANCE[:,0], "%s_neighbors_distance_0.h5"%nt.this_looper.out_prefix, 'distance_0')

        dump_core_vals(cores,np.log10(nt.mean_density), "%s_neighbors_mean_density_log.h5"%nt.this_looper.out_prefix,'mean_density')
        dump_core_vals(cores,np.log10(nt.NEXT_DISTANCE[:,0]), "%s_neighbors_distance_0_log.h5"%nt.this_looper.out_prefix, 'distance_0')

        tt=densities[i]

        G=1620./(4*np.pi)
        tff_global = np.sqrt(3*np.pi/(32*G*1))

        outname = tt.this_looper.out_prefix
        dump_core_vals( tt.cores_used,tt.collapse_times, fname='plots_to_sort/%s_first_collapse.h5'%outname,setname='collapse_times')
        dump_core_vals( tt.cores_used,np.log10(tt.tc_mean), fname='plots_to_sort/%s_mean_collapse_log.h5'%outname,setname='collapse_times')
        dump_core_vals( tt.cores_used,np.log10(tt.collapse_times), fname='plots_to_sort/%s_ct_log.h5'%outname,setname='collapse_times')
        dump_core_vals( tt.cores_used,np.log10(tt.collapse_times/tff_global), fname='plots_to_sort/%s_ct_glob_log.h5'%outname,setname='collapse_times')
        dump_core_vals( tt.cores_used,np.log10(nar(tt.collapse_times)/nar(tt.tff_local_list)), 
                       fname='plots_to_sort/%s_ct_local_log.h5'%outname,setname='collapse_times')


