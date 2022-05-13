
from starter2 import *

import three_loopers_u500 as TL
import movie_frames 
import heat_map

def simple_v(this_looper,core_list=None, symlog=False):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    mask = movie_frames.quantized_mask(this_looper).flatten()
    times=thtr.times[mask]+0 #the zero makes a copy
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

        vx = ms.raw_vx[:,mask]
        vy = ms.raw_vx[:,mask]
        vz = ms.raw_vx[:,mask]
        vr = ms.rel_vmag[:,mask]
        vr_extents(vr)
        vvv = [vx,vy,vz,vr]
        Nlin = 32
        Nlog = 32
        if symlog:
            bins1 = np.geomspace(1,16,Nlog)
            bins2 = np.linspace(-1,1,Nlin)
            bins = np.unique(np.concatenate( [ -bins1[::-1], bins2, bins1]))
        else:
            bins = np.linspace(-16,16,65)
        #fig2,ax2=plt.subplots(1,1)
        #ax2.plot(bins)
        #fig2.savefig('/home/dccollins/PigPen/bins.png')

        for nv,v in enumerate(vvv):
            label = [r'$v_x$',r'$v_y$',r'$v_z$',r'$v_r$'][nv]



            stuff = heat_map.heat_map( v, times, ax=axlist[nv],bins=bins)
            axlist[nv].plot(times, times*0, c='k')
            axlist[nv].plot(times, times*0+1, c='k')
            axlist[nv].plot(times, times*0-1, c='k')
            axbonk(axlist[nv],xlabel=r'$t/t_{ff}$', ylabel=label)
            if symlog:
                axlist[nv].set_yscale('symlog',linthresh=1)



        if symlog:
            logornot="_symlog"
        else:
            logornot=""
        outname='plots_to_sort/%s_v_t_heat_c%04d%s.png'%(this_looper.sim_name,core_id,logornot)
        fig.savefig(outname)
        print(outname)
    print("VR extents",vr_extents.minmax)

sims=['u501', 'u502','u503']
for sim in sims:
    simple_v(TL.loops[sim],symlog=False)#, core_list=[149,323])
    #break

