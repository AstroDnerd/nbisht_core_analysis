from starter2 import *
import matplotlib.image as mpimg

from scipy.optimize import curve_fit
import data_locations as dl
reload(dl)
plt.close('all')


class trial():
    def __init__(self):
        self.collapse_times=[]
        self.used_cores=[]
        self.tff_local_list=[]
        self.rho_0_list = []
        self.rho_c_list = []
        self.tc=[]
        self.tc_mean=[]
        self.rho_col_mean=[]
        self.rho_col_first=[]
    def count_particles(self,the_looper):
        cores=[]
        nparticles=[]
        for core in the_looper.target_indices:
            cores.append(core)
            nparticles.append(the_looper.target_indices[core].size)
        return {'core':nar(core),'nparticles':nar(nparticles)}
    def run(self,this_looper,do_all_plots=True,core_list=None):

        thtr = this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        for core_id in core_list:
            print('go ', core_id)
            self.used_cores.append(core_id)
            n0=0

            fig,ax=plt.subplots(1,1)

            ms = trackage.mini_scrubber(thtr,core_id)
            density = thtr.c([core_id],'density')
            cell_volume = thtr.c([core_id],'cell_volume')
            rho_mean = density.mean(axis=0)[:]

            #collapse time
            if (rho_mean>1e3).any():
                first_collapse_index = np.where( (density>1e3).sum(axis=0) >0 )[0][0] 

                t_collapse = thtr.times[first_collapse_index]
                rho_col = density[:,first_collapse_index].max()

                mean_collapse_index= np.where( rho_mean>1e3 )[0][0] 
                tc_mean=thtr.times[mean_collapse_index]
                rho_col_mean=density[:,mean_collapse_index].max()

                if do_all_plots:
                    ax.plot( [tc_mean-0.005, tc_mean+0.005],[rho_col_mean,rho_col_mean],c=[0.5]*3)
                    ax.plot( [tc_mean      , tc_mean     ], [rho_col_mean/5,5*rho_col_mean],c=[0.5]*3)
                    print('WWWW ',tc_mean, rho_col_mean)


            else:
                t_collapse = -1
                rho_col = -1
                tc_mean = -1
                rho_col_mean = -1

            self.collapse_times.append(t_collapse)
            self.tc_mean.append(tc_mean)
            self.rho_col_mean.append(rho_col_mean)
            self.rho_col_first.append(rho_col)

            t0 = thtr.times[0]
            t1 = thtr.times[-1]
            rho0 = (density[:,0]*cell_volume[:,0]).sum()/cell_volume[:,0].sum()
            rho1 = density[:,0].max() 
            alpha = 1.8
            G=1620./(4*np.pi)

            tff_global = np.sqrt(3*np.pi/(32*G*1))
            tff_local = np.sqrt(3*np.pi/(32*G*rho0))

            self.tff_local_list.append(tff_local)
            self.rho_0_list.append(rho0)

            if do_all_plots:
                norm = mpl.colors.Normalize()
                norm.autoscale( np.log10(density[:,n0]))
                cmap = mpl.cm.jet
                color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)
                if ms.nparticles<100:
                    this_slice=slice(None)
                else:
                    this_slice=slice(None,None,10)

                for npart in list(range(ms.nparticles))[this_slice]:
                    c = color_map.to_rgba(np.log10(density[npart,n0]))
                    ax.plot( thtr.times, density[npart,:],c=c,linewidth=.1)#linestyle=':')
                ax.plot(tsorted, rho_mean,c='k')

                if rho_col > 0:
                    tc =t_collapse*(1-(rho_col/rho0)**(-1./alpha))**-0.5
                    rho_c = 3*np.pi/(32*G*tc**2)
                    rho_tc = rho_c*(1-(tsorted/tc)**2)**-alpha
                    ax.plot( tsorted, rho_tc, c='g')
                    print("stuff rhoc/rho0 %0.2f"%(rho_c/rho0))
                    rho_tff = rho0*(1-(tsorted/tff_local)**2)**-alpha
                    ax.plot( tsorted, rho_tff, c='b')
                    #rho_c = 3*np.pi/(32*G*tc**2)
                    self.tc.append(tc)
                    self.rho_c_list.append(rho_c)
                else:
                    self.tc.append(-1)
                    self.rho_c_list.append(-1)

                ax.legend(loc=0)
                ylim = [thtr.track_dict['density'].min(), thtr.track_dict['density'].max()]
                xlim = [0, 0.05]

                axbonk(ax,xscale='linear',yscale='log',xlabel='t',ylabel=r'$\rho$', ylim=ylim,xlim=xlim)
                oname = "%s/%s_density_6_c%04d"%(dl.output_directory,this_looper.out_prefix,core_id)
                fig.savefig(oname)
                print("Saved "+oname)

out_prefix="u05u10u11"

if 'looper1' not in dir():
    file_list=glob.glob(dl.every_ten['u05'])
    looper1=looper.core_looper(directory=dl.sims['u05'])
    for nfile,fname in enumerate(file_list):
        looper1.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    thtr = looper1.tr
    looper1.out_prefix='u05'
    thtr.sort_time()

if 'looper2' not in dir():
    looper2=looper.core_looper(directory=dl.sims['u10'])
    file_list=glob.glob(dl.every_ten['u10'])
    for nfile,fname in enumerate(file_list):
        looper2.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    looper2.out_prefix='u10'

if 'looper3' not in dir():
    looper3=looper.core_looper(directory=dl.sims['u11'])
    file_list=glob.glob(dl.every_ten['u11'])
    for nfile,fname in enumerate(file_list):
        looper3.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    looper3.out_prefix='u11'

if 'tool1' not in dir():
    tool1=trial()
    tool1.run(looper1,do_all_plots=True)
    tool2=trial()
    tool2.run(looper2,do_all_plots=True)
    tool3=trial()
    tool3.run(looper3,do_all_plots=True)


if 0:
    def dump_core_vals(cores,vals,fname='out.h5',setname='values'):
        fptr=h5py.File(fname,'w')
        fptr.create_dataset("core_ids",data=cores)
        fptr.create_dataset(setname,data=vals)
        fptr.close()

    #dump_core_vals( tool2.core_list, tool2.collapse_times, fname='plots_to_sort/u05_ct.h5',setname='collapse_times')
    #dump_core_vals( tool2.core_list, tool2.collapse_times/tff_global, fname='plots_to_sort/u05_ct_glob.h5',setname='collapse_times')
    #dump_core_vals( tool2.core_list, nar(tool2.collapse_times)/nar(tool2.tff_local_list), 
    #               fname='plots_to_sort/u05_ct_local.h5',setname='collapse_times')


#stuff1 = tool1.count_particles( looper1)
#stuff2 = tool2.count_particles( looper2)
#stuff3 = tool3.count_particles( looper3)

if 0:
    G=1620./(4*np.pi)
    tff_global = np.sqrt(3*np.pi/(32*G*1))
    fig, ax = plt.subplots(1,1) 
    ct1=nar(tool1.collapse_times)
    ct2=nar(tool2.collapse_times)
    ct3=nar(tool3.collapse_times)
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

if 1:
    G=1620./(4*np.pi)
    tff_global = np.sqrt(3*np.pi/(32*G*1))
    fig, ax = plt.subplots(1,1) 
    Lct1=nar(tool1.tc_mean)
    Lct2=nar(tool2.tc_mean)
    Lct3=nar(tool3.tc_mean)
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

    ax.scatter(stuff1['nparticles'][ok1],ct1[ok1]/nar(tool1.tff_local_list)[ok1],label='u05')
    ax.scatter(stuff2['nparticles'][ok2],ct2[ok2]/nar(tool2.tff_local_list)[ok2],label='u10')
    ax.scatter(stuff3['nparticles'][ok3],ct3[ok3]/nar(tool3.tff_local_list)[ok3],label='u11')
    axbonk(ax,xlabel='N',ylabel=r'$t_c/t_{\rm{ff,local}}$',xscale='log')
    ax.legend(loc=0)
    oname = "%s/%s_np_tc_local"%(dl.output_directory,out_prefix)
    fig.savefig( oname)
    plt.close(fig)
    print(oname)

if 0:

    fig, ax = plt.subplots(1,1) 

    ax.scatter(stuff1['nparticles'][ok1],ct1[ok1]/nar(tool1.tff_harm_list)[ok1],label='u05')
    ax.scatter(stuff2['nparticles'][ok2],ct2[ok2]/nar(tool2.tff_harm_list)[ok2],label='u10')
    ax.scatter(stuff3['nparticles'][ok3],ct3[ok3]/nar(tool3.tff_harm_list)[ok3],label='u11')
    axbonk(ax,xlabel='N',ylabel=r'$t_c/t_{\rm{ff,harmonic}}$',xscale='log')
    ax.legend(loc=0)
    oname = "%s/%s_np_tc_harmonic"%(dl.output_directory,out_prefix)
    fig.savefig( oname)
    plt.close(fig)
    print(oname)

if 0:

    fig, ax = plt.subplots(1,1) 

    rho_a =min(tool1.rho_0_list)
    rho_b =max(tool1.rho_0_list)
    ax.plot( [rho_a,rho_b],[rho_a,rho_b], c='k')
    ax.plot( [2*rho_a,2*rho_b],[rho_a,rho_b], c='k')
    ax.scatter(tool1.rho_0_list, tool1.rho_harmonic,label='u05',marker='.',s=.9,c='r')
    ax.scatter(tool2.rho_0_list, tool2.rho_harmonic,label='u10',marker='.',s=.9,c='g')
    ax.scatter(tool3.rho_0_list, tool3.rho_harmonic,label='u11',marker='.',s=.9,c='b')
    axbonk(ax,xlabel=r'$\rho_0$',ylabel=r'$\rho_{\rm{harm}}$',xscale='log',yscale='log')
    ax.legend(loc=0)
    oname = "%s/%s_rho_harmonic"%(dl.output_directory,out_prefix)
    fig.savefig( oname)
    
    print(oname)
