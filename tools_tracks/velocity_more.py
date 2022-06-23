
from starter2 import *

import three_loopers_u500 as TL
import movie_frames 

class vstuff():
    def __init__(self,this_looper):
        self.this_looper=this_looper
    
    def run(self,core_list=None, symlog=False, do_plots=False, atool=None):
        this_looper=self.this_looper

        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        thtr=this_looper.tr
        mask = movie_frames.quantized_mask(this_looper).flatten()
        times=thtr.times[mask]+0 #the zero makes a copy
        times.shape=times.size,1
        times=times/colors.tff
        vr_extents=extents()

        if do_plots:
            fig,ax=plt.subplots(2,2,figsize=(12,12))
            axlist=ax.flatten()
            fig.subplots_adjust(wspace=0, hspace=0)
            ax[0][0].xaxis.tick_top()
            ax[0][1].xaxis.tick_top()
            ax[0][1].xaxis.set_label_position('top')
            ax[0][0].xaxis.set_label_position('top')
            ax[1][1].yaxis.tick_right()
            ax[0][1].yaxis.tick_right()
            ax[1][1].yaxis.set_label_position('right')
            ax[0][1].yaxis.set_label_position('right')

            ax5=ax[0][0].twinx()

        self.sigma_3d=np.zeros([len(core_list),len(times)])

        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True)

            if do_plots:
                for aaa in ax.flatten():
                    aaa.clear()
                ax5.clear()
            if ms.nparticles < 1000:
                sl=slice(None)
                c=[0.5]*4
            else:
                sl = slice(None,None,10)
                #c=[0,0,0,0.1]
                c=[0.1]*4

            #ms.particle_pos(core_id)
            ms.get_central_at_once(core_id)

            vx = ms.rel_vx[sl].transpose()[mask,:]
            vy = ms.rel_vy[sl].transpose()[mask,:]
            vz = ms.rel_vz[sl].transpose()[mask,:]
            vr = ms.rel_vmag[sl].transpose()[mask,:]
            vr_extents(vr)
            vvv = [vx,vy,vz,vr]
            if do_plots:
                for nv,v in enumerate(vvv):
                    axlist[nv].plot(times , v, c=c, linewidth=0.1)
                    label = [r'$v_x$',r'$v_y$',r'$v_z$',r'$|(v-\bar{v})|$'][nv]
                    vext = [-20,20]

                    axlist[nv].plot(times, times*0, c='k',linewidth=0.2)
                    axlist[nv].plot(times, times*0+1, c='k',linewidth=0.2)
                    axlist[nv].plot(times, times*0-1, c='k',linewidth=0.2)
                    axbonk(axlist[nv],xlabel=r'$t/t_{ff}$', ylabel=label, ylim=vext)
                    if symlog:
                        axlist[nv].set_yscale('symlog',linthresh=1)
                axlist[3].set_ylim([0,vext[1]])

                if atool is not None:
                    ax5.plot( times, atool.Alpha_rho_r[nc,:], c='r',label=r'$\alpha$')
                    ax5.plot( times, atool.Pearson_rho_r[nc,:], c='g',label=r'$r-\rho$')
                    ax5.plot( times, atool.Pearson_GE_r[nc,:], c='b',label=r'$r-GE$')
                    ax5.set_ylim(-2.5,2.5)
                    ax5.legend(loc=2)


                logornot=""
                if symlog:
                    logornot="_symlog"



                outname='plots_to_sort/%s_v_t_c%04d.png'%(this_looper.sim_name,core_id)
                fig.savefig(outname)
                print(outname)

sims=['u501', 'u502','u503']
import atool
reload(atool)
if 'alpha_time' not in dir():
    alpha_time={}
    for sim in sims:
        alpha_time[sim] = atool.atool(TL.loops[sim])
        alpha_time[sim].run()

sims=['u501', 'u502','u503']
for sim in sims:
    vs = vstuff(TL.loops[sim])
    vs.run(symlog=False, do_plots=True, atool=alpha_time[sim])#, core_list=[0,87,149,323])

