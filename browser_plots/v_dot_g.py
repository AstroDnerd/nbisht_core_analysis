
from starter2 import *

import three_loopers_u500 as TL
import movie_frames 

def simple_v(this_looper,core_list=None, symlog=False):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    mask = movie_frames.quantized_mask(this_looper).flatten()
    times=thtr.times[mask]+0 #the zero makes a copy
    times.shape=times.size,1
    times=times/colors.tff
    vr_extents=extents()
    fig,ax=plt.subplots(2,2,figsize=(12,12))
    fig.subplots_adjust(wspace=0, hspace=0)
    ax[0][0].xaxis.tick_top()
    ax[0][1].xaxis.tick_top()
    ax[0][1].xaxis.set_label_position('top')
    ax[0][0].xaxis.set_label_position('top')
    ax[1][1].yaxis.tick_right()
    ax[0][1].yaxis.tick_right()
    ax[1][1].yaxis.set_label_position('right')
    ax[0][1].yaxis.set_label_position('right')
    for core_id in core_list:
        axlist=ax.flatten()
        for aaa in axlist:
            aaa.clear()

            
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True)
        ms.particle_pos(core_id)

        if ms.nparticles < 1000:
            sl=slice(None)
            c=[0.5]*4
        else:
            sl = slice(None,None,10)
            #c=[0,0,0,0.1]
            c=[0.1]*4

        vx = ms.raw_vx[sl].transpose()[mask,:]
        vy = ms.raw_vy[sl].transpose()[mask,:]
        vz = ms.raw_vz[sl].transpose()[mask,:]
        vr = ms.rel_vmag[sl].transpose()[mask,:]
        vr_extents(vr)
        vvv = [vx,vy,vz,vr]
        fig2,ax2=plt.subplots(1,1)
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

            ax2.hist(v[-1,:],histtype='step',label=label)
        ax2.legend(loc=0)
        logornot=""
        if symlog:
            logornot="_symlog"
        fig2.savefig('plots_to_sort/vhist_%s%s.png'%(this_looper.sim_name,logornot))



        outname='plots_to_sort/%s_v_t_c%04d.png'%(this_looper.sim_name,core_id)
        fig.savefig(outname)
        print(outname)
    print("VR extents",vr_extents.minmax)

sims=['u501', 'u502','u503']
for sim in sims:
    simple_v(TL.loops[sim],symlog=False)#, core_list=[149,323])
    #break

