from starter2 import *
import matplotlib.image as mpimg

from scipy.optimize import curve_fit

import density_tools
reload(density_tools)

import data_locations as dl
reload(dl)
plt.close('all')


import three_loopers_mountain_top as TLM

if 1:
    dt3={}
    for this_simname in ['u301','u302','u303']:
        dt3[this_simname] = density_tools.trial( TLM.loops[this_simname])
        dt3[this_simname].run()



if 'do_all_plots' not in dir():
    do_all_plots = False


if 'dt14' not in dir() and False:
    print('wut')
    import looper_u14
    reload(looper_u14)

    dt14=trial(looper_u14.looper14)
    dt14.run(do_all_plots=True)

if do_all_plots and False:
    out_prefix="u05u10u11"

    import three_loopers as tl

    if 'density_tool1' not in dir():
        density_tool1=trial(tl.looper1)
        density_tool1.run(do_all_plots=False)
        density_tool2=trial(tl.looper2)
        density_tool2.run(do_all_plots=False)
        density_tool3=trial(tl.looper3)
        density_tool3.run(do_all_plots=False)


    if 0:
        def dump_core_vals(cores,vals,fname='out.h5',setname='values'):
            fptr=h5py.File(fname,'w')
            fptr.create_dataset("core_ids",data=cores)
            fptr.create_dataset(setname,data=vals)
            fptr.close()

        for tt, ll in [ [density_tool1,looper1],[density_tool2,looper2],[density_tool3,looper3]]:
            outname=ll.out_prefix
            dump_core_vals( tt.core_list, tt.collapse_times, fname='plots_to_sort/%s_ct.h5'%outname,setname='collapse_times')
            dump_core_vals( tt.core_list, tt.collapse_times/tff_global, fname='plots_to_sort/%s_ct_glob.h5'%outname,setname='collapse_times')
            dump_core_vals( tt.core_list, nar(tt.collapse_times)/nar(tt.tff_local_list), 
                           fname='plots_to_sort/%s_ct_local.h5'%outname,setname='collapse_times')


#stuff1 = density_tool1.count_particles( looper1)
#stuff2 = density_tool2.count_particles( looper2)
#stuff3 = density_tool3.count_particles( looper3)

    if 0:
        G=1620./(4*np.pi)
        tff_global = np.sqrt(3*np.pi/(32*G*1))
        fig, ax = plt.subplots(1,1) 
        ct1=nar(density_tool1.collapse_times)
        ct2=nar(density_tool2.collapse_times)
        ct3=nar(density_tool3.collapse_times)
        ok1=ct1>0
        ok2=ct2>0
        ok3=ct3>0
        ax.hist( ct1[ok1], histtype='step',label='u05')
        ax.hist( ct2[ok2], histtype='step',label='u10')
        ax.hist( ct3[ok3], histtype='step',label='u11')
        axbonk(ax, xlabel=r'$t_c/t_{\rm{ff}}$', ylabel = 'N')
        oname = "%s/%s_tc_hist"%(dl.output_directory,out_prefix)
        fig.savefig( oname)
        plt.close(fig)
        print(oname)

    if 0:
        G=1620./(4*np.pi)
        tff_global = np.sqrt(3*np.pi/(32*G*1))
        fig, ax = plt.subplots(1,1) 
        Lct1=nar(density_tool1.tc_mean)
        Lct2=nar(density_tool2.tc_mean)
        Lct3=nar(density_tool3.tc_mean)
        Lok1=Lct1>0
        Lok2=Lct2>0
        Lok3=Lct3>0
        ax.hist(Lct1[Lok1]/tff_global, histtype='step',label='u05',color='r')
        ax.hist(Lct2[Lok2]/tff_global, histtype='step',label='u10',color='g')
        ax.hist(Lct3[Lok3]/tff_global, histtype='step',label='u11',color='b')
        ax.scatter(looper1.tr.times.max()/tff_global,20,color='r')
        ax.scatter(looper2.tr.times.max()/tff_global,20,color='g')
        ax.scatter(looper3.tr.times.max()/tff_global,20,color='b')
        axbonk(ax, xlabel=r'$t_{c,mean}/t_{\rm{ff}}$', ylabel = 'N')
        oname = "%s/%s_tc_mean_hist"%(dl.output_directory,out_prefix)
        fig.savefig( oname)
        plt.close(fig)
        print(oname)

    if 0:
        fig, ax = plt.subplots(1,1) 
        ax.scatter(stuff1['nparticles'][ok1],ct1[ok1]/tff_global,label='u05')
        ax.scatter(stuff2['nparticles'][ok2],ct2[ok2]/tff_global,label='u10')
        ax.scatter(stuff3['nparticles'][ok3],ct3[ok3]/tff_global,label='u11')
        axbonk(ax,xlabel='N',ylabel=r'$t_c/t_{\rm{ff}}$',xscale='log')
        ax.legend(loc=0)
        oname = "%s/%s_np_tc_hist"%(dl.output_directory,out_prefix)
        fig.savefig( oname)
        plt.close(fig)
        print(oname)

    if 0:

        fig, ax = plt.subplots(1,1) 

        ax.scatter(stuff1['nparticles'][ok1],ct1[ok1]/nar(density_tool1.tff_local_list)[ok1],label='u05')
        ax.scatter(stuff2['nparticles'][ok2],ct2[ok2]/nar(density_tool2.tff_local_list)[ok2],label='u10')
        ax.scatter(stuff3['nparticles'][ok3],ct3[ok3]/nar(density_tool3.tff_local_list)[ok3],label='u11')
        axbonk(ax,xlabel='N',ylabel=r'$t_c/t_{\rm{ff,local}}$',xscale='log')
        ax.legend(loc=0)
        oname = "%s/%s_np_tc_local"%(dl.output_directory,out_prefix)
        fig.savefig( oname)
        plt.close(fig)
        print(oname)

    if 0:

        fig, ax = plt.subplots(1,1) 

        ax.scatter(stuff1['nparticles'][ok1],ct1[ok1]/nar(density_tool1.tff_harm_list)[ok1],label='u05')
        ax.scatter(stuff2['nparticles'][ok2],ct2[ok2]/nar(density_tool2.tff_harm_list)[ok2],label='u10')
        ax.scatter(stuff3['nparticles'][ok3],ct3[ok3]/nar(density_tool3.tff_harm_list)[ok3],label='u11')
        axbonk(ax,xlabel='N',ylabel=r'$t_c/t_{\rm{ff,harmonic}}$',xscale='log')
        ax.legend(loc=0)
        oname = "%s/%s_np_tc_harmonic"%(dl.output_directory,out_prefix)
        fig.savefig( oname)
        plt.close(fig)
        print(oname)

    if 0:

        fig, ax = plt.subplots(1,1) 

        rho_a =min(density_tool1.rho_0_list)
        rho_b =max(density_tool1.rho_0_list)
        ax.plot( [rho_a,rho_b],[rho_a,rho_b], c='k')
        ax.plot( [2*rho_a,2*rho_b],[rho_a,rho_b], c='k')
        ax.scatter(density_tool1.rho_0_list, density_tool1.rho_harmonic,label='u05',marker='.',s=.9,c='r')
        ax.scatter(density_tool2.rho_0_list, density_tool2.rho_harmonic,label='u10',marker='.',s=.9,c='g')
        ax.scatter(density_tool3.rho_0_list, density_tool3.rho_harmonic,label='u11',marker='.',s=.9,c='b')
        axbonk(ax,xlabel=r'$\rho_0$',ylabel=r'$\rho_{\rm{harm}}$',xscale='log',yscale='log')
        ax.legend(loc=0)
        oname = "%s/%s_rho_harmonic"%(dl.output_directory,out_prefix)
        fig.savefig( oname)
        
        print(oname)
