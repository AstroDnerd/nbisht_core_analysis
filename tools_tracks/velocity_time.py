from starter2 import *
import data_locations as dl
import davetools
reload(looper)
reload(trackage)
reload(davetools)

plt.close('all')

class tool_velocity():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.rho_extents = None
        self.r_extents = None
        self.vel_ext= None

    def run(self,core_list=None):
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()
        rm = rainbow_map(len(all_cores))
        if self.rho_extents is None:
            self.rho_extents=davetools.extents()
            self.r_extents=davetools.extents()
            for nc,core_id in enumerate(all_cores):
                ms = trackage.mini_scrubber(thtr,core_id)
                if ms.nparticles == 1:
                    continue
                density = thtr.c([core_id],'density')
                self.rho_extents(density)
                self.r_extents(ms.r)

        plt.close('all')
        fig, axd1=plt.subplots(1,1)
        fig4, ax4=plt.subplots(2,2)
        ax40 = ax4[0][0]
        ax41 = ax4[0][1]
        ax42 = ax4[1][0]
        ax43 = ax4[1][1]

        if self.vel_ext is None:
            #If we haven't done it yet, we'll need to run the velocity extents.
            self.vel_ext = extents()
            print("doing velocity extents")
            for nc,core_id in enumerate(core_list):
                ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
                sigma2 = np.sqrt((ms.rel_vx**2 + ms.rel_vy**2 + ms.rel_vx**2).mean(axis=0))
                self.vel_ext(sigma2)
                self.vel_ext(ms.rel_vx)
                self.vel_ext(ms.rel_vy)
                self.vel_ext(ms.rel_vz)
                self.vel_ext(ms.rel_vmag)

        print('make plots')
        for nc,core_id in enumerate(core_list):

            #miniscrubber computes distance, r^2, several other quantities
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)

    
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

            if ms.nparticles<100:
                this_slice=slice(None)
            else:
                this_slice=slice(None,None,10)
            for npart in list(range(ms.nparticles))[this_slice]:
                ax40.set_title("%s particles"%ms.nparticles)
                this_r=ms.r[npart,:]+0
                this_r[ this_r < 1./2048] = 1./2048
                r_un = nar(sorted(np.unique(this_r)))
                this_t = thtr.times
                c=[0.5]*3
                outname4 = '%s/%s_vi_t_rel_volume_centroid_c%04d'%(dl.output_directory,self.this_looper.out_prefix,core_id)
                ax40.plot(this_t, ms.rel_vx[npart,:], c=c,lw=.1)
                ax41.plot(this_t, ms.rel_vy[npart,:], c=c,lw=.1)
                ax42.plot(this_t, ms.rel_vz[npart,:], c=c,lw=.1)
                
                ax43.plot(this_t, ms.rel_vmag[npart,:], c=c,lw=.1)


            mean_vmag = np.mean(ms.rel_vmag,axis=0)
            sigma2 = np.sqrt((ms.rel_vx**2 + ms.rel_vy**2 + ms.rel_vx**2).mean(axis=0))
            ax43.plot(this_t, mean_vmag, c='r')
            ax43.plot(this_t, sigma2, 'r:')

            unique_mask = ms.compute_unique_mask(core_id, dx=1./2048)
            um = unique_mask
            mean_vmag = np.mean(ms.rel_vmag[um],axis=0)
            sigma2 = np.sqrt((ms.rel_vx[um]**2 + ms.rel_vy[um]**2 + ms.rel_vx[um]**2).mean(axis=0))
            ax43.plot(this_t, mean_vmag, c='g')
            ax43.plot(this_t, sigma2, 'g:')

            ax43.plot(this_t, np.mean(np.abs(ms.vr_rel),axis=0), c='b')

            #axd1.plot([ms.r.min(),ms.r.max()], [1,1], c=[0.5]*4)
            limits = np.abs(self.vel_ext.minmax).max()
            limits = [-limits,limits]

            labs = ['vx','vy','vz','vtotal']
            for iii,ax4i in enumerate([ax40,ax41,ax42,ax43]):
                davetools.axbonk(ax4i,xlabel='r',ylabel=labs[iii],xscale='linear',ylim=limits)
                ax4i.set_yscale('symlog',linthresh=2)
                #ax4i.set_yscale('linear')
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
                axd1.set_xscale('symlog',linthresh=1./2048)
                axd1.set_xlim([0,1])
                axd1.set_yscale('symlog',linthresh=1.)
                axd1.set_ylim([0,15])
                outname = '%s/vr_vtotal_c%04d'%(dl.output_directory,core_id)

            #outname = '%s/vr_vtotal_c%04d'%(dl.output_directory,core_id)
            #fig.savefig(outname)
            #print("saved "+outname)
            fig4.savefig(outname4)
            print("saved "+outname4)
            plt.close(fig4)

import three_loopers_1tff as TL

tv1 = tool_velocity( TL.looper1 )
#tv2 = tool_velocity( TL.looper2 )
#tv3 = tool_velocity( TL.looper3 )
tv1.run(core_list=[32])
#tv2.run()
#tv3.run()
