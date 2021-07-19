
from starter2 import *
import matplotlib.image as mpimg

from scipy.optimize import curve_fit
import data_locations as dl
reload(dl)
plt.close('all')


class trial():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.collapse_times=[]
        self.cores_used=[]
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
    def run(self,do_all_plots=True,core_list=None):

        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            if ms.nparticles < 3:
                continue
            print('go ', core_id)
            self.cores_used.append(core_id)
            n0=0

            fig,ax=plt.subplots(1,1)

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
                xlim = [0,thtr.times.max()]

                axbonk(ax,xscale='linear',yscale='log',xlabel='t',ylabel=r'$\rho$', ylim=ylim,xlim=xlim)
                oname = "%s/%s_density_6_c%04d"%(dl.output_directory,self.this_looper.out_prefix,core_id)
                print(oname)
                fig.savefig(oname)
                print("Saved "+oname)
