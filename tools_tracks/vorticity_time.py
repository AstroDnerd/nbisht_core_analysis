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
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            if ms.nparticles < 3:
                continue
            print('go ', core_id)
            self.used_cores.append(core_id)
            n0=0

            fig,axes=plt.subplots(2,2, sharex=True)
            #fig.subplots_adjust(wspace=0, hspace=0)
            ax1 = axes[0][0]
            ax2 = axes[1][0]
            ax3 = axes[0][1]
            ax4 = axes[1][1]

            vorticity = thtr.c([core_id],'vorticity_magnitude')
            density = thtr.c([core_id],'density')
            cell_volume = thtr.c([core_id],'cell_volume')
            Omegax = ((ms.this_y*ms.rel_vz*-ms.this_z*ms.rel_vy)*cell_volume).sum(axis=0)/cell_volume.sum()/ms.r2
            Omegay = ((ms.this_z*ms.rel_vx*-ms.this_x*ms.rel_vz)*cell_volume).sum(axis=0)/cell_volume.sum()/ms.r2
            Omegaz = ((ms.this_x*ms.rel_vy*-ms.this_y*ms.rel_vx)*cell_volume).sum(axis=0)/cell_volume.sum()/ms.r2
            Omegamag = np.sqrt(Omegax**2+Omegay**2+Omegaz**2)

            I = (density*cell_volume*ms.r**2).sum(axis=0)

            plt.clf()
            plt.scatter(density[:,12], Omegamag[:,12]**2)
            plt.yscale('log')
            plt.xscale('log')
            plt.savefig('/home/dccollins/PigPen/vort_test.png')

            jx = cell_volume*density*(ms.this_y*ms.rel_vz*-ms.this_z*ms.rel_vy)
            jy = cell_volume*density*(ms.this_z*ms.rel_vx*-ms.this_x*ms.rel_vz)
            jz = cell_volume*density*(ms.this_x*ms.rel_vy*-ms.this_y*ms.rel_vx)
            j2 = jx*jx+jy*jy+jz*jz
            jmag = np.sqrt(j2)
            #rho_mean = density.mean(axis=0)[:]

            ext_j = extents()
            if do_all_plots:
                j_scale = 10**np.mean(np.log10(jmag[:,0]))
                const = [j_scale]*len(thtr.times)
                ax1.plot(thtr.times, const,c=[0.7]*4)
                omega_scale = 10**np.mean(np.log10(vorticity[:,0]))
                const = [omega_scale]*len(thtr.times)
                ax2.plot(thtr.times, const,c=[0.7]*4)
                #norm = mpl.colors.Normalize()
                #norm.autoscale( np.log10(density[:,n0]))
                #cmap = mpl.cm.jet
                #color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)
                ax3.plot( thtr.times, Omegamag,c='k',linewidth=.1)#linestyle=':')
                ax4.plot( thtr.times, I,c='k',linewidth=.1)#linestyle=':')
                if ms.nparticles<100:
                    this_slice=slice(None)
                else:
                    this_slice=slice(None,None,10)

                for npart in list(range(ms.nparticles))[this_slice]:
                    #c = color_map.to_rgba(np.log10(density[npart,n0]))
                    c=[0.5]*4
                    #ax.plot( thtr.times, vorticity[npart,:],c=c,linewidth=.1)#linestyle=':')
                    this_j = jmag[npart,:]
                    rel_j = np.abs(1 - this_j/j_scale)
                    if np.isnan(rel_j).any() or np.isinf(rel_j.any()):
                        pdb.set_trace()
                    ext_j(this_j)
                    ax1.plot( thtr.times, this_j ,c=c,linewidth=.1)#linestyle=':')
                    ax2.plot( thtr.times, vorticity[npart,:],c=c,linewidth=.1)#linestyle=':')

                for nquan in [0,1]:
                    this_ax = [ax1,ax2][nquan]
                    val = [jmag, vorticity][nquan]
                    this_ax.plot( thtr.times, val.mean(axis=0),c='k')
                    p25 = np.percentile(val,25,axis=0)
                    this_ax.plot( thtr.times, p25,c=[0.5,0,0,0.5])
                    p75 = np.percentile(val,25,axis=0)
                    this_ax.plot( thtr.times, p75,c=[0.5,0,0,0.5])


                ax1.plot( thtr.times, this_j ,c=c,linewidth=.1)#linestyle=':')
                ax2.plot( thtr.times, vorticity[npart,:],c=c,linewidth=.1)#linestyle=':')
                #ax.plot(tsorted, rho_mean,c='k')


                #ylim1 = [np.sqrt(j2).min(), np.sqrt(j2).max()]
                ylim1 = ext_j.minmax
                ylim2 = [vorticity.min(), vorticity.max()]
                xlim = [0, 0.05]

                axbonk(ax1,xscale='linear',yscale='log',xlabel='t',ylabel=r'$j=\rho r \times v$', ylim=ylim1,xlim=xlim)
                axbonk(ax2,xscale='linear',yscale='log',xlabel='t',ylabel=r'$\omega=\nabla \times v$', ylim=ylim2,xlim=xlim)
                axbonk(ax3,xscale='linear',yscale='log',xlabel='t',ylabel=r'$\Omega = r \times v$')
                axbonk(ax4,xscale='linear',yscale='log',xlabel='t',ylabel=r'$I =\Sigma m r^2$')
                oname = "%s/%s_vorticity_c%04d"%(dl.output_directory,this_looper.out_prefix,core_id)
                fig.savefig(oname)
                print("Saved "+oname)

out_prefix="u05u10u11"
import three_loopers_energy as tle  
reload(tle)
if 'clobber' not in dir():
    clobber=False
if 'tool1' not in dir() or clobber:
    tool1=trial()
    tool1.run(tle.looper1,do_all_plots=True,core_list=[14])
#   tool2=trial()
#   tool2.run(tle.looper2,do_all_plots=True)
#   tool3=trial()
#   tool3.run(tle.looper3,do_all_plots=True)


