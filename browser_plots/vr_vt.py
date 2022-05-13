
from starter2 import *

import three_loopers_u500 as TL
import movie_frames 
import heat_map

def vr_vt(this_looper,core_list=None, symlog=False):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    mask = movie_frames.quantized_mask(this_looper).flatten()
    times=thtr.times[mask]+0 #the zero makes a copy
    times=times/colors.tff
    vr_extents=extents()
    fig,ax=plt.subplots(1,2,figsize=(12,8))
    axlist=ax.flatten()

    for core_id in core_list:
        for aaa in axlist:
            aaa.clear()
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True)

        vr = ms.vr_rel[:,mask]
        vt = np.sqrt(ms.vt2_rel[:,mask])
        #r  = ms.r[:,mask]
        #omega = vt/r
        #print(omega.min(),omega.max())
        Nlin = 32
        Nlog = 32
        if symlog:
            bins1 = np.geomspace(1,20,Nlog)
            bins2 = np.linspace(-1,1,Nlin)
            bins = np.unique(np.concatenate( [ -bins1[::-1], bins2, bins1]))
        else:
            bins=np.linspace(-16,16,65)

        vvv=[vr,vt]
        for nv,v in enumerate(vvv):
            label = [r'$v_r$',r'$v_t$',r'$\omega$'][nv]
            stuff = heat_map.heat_map( v, times, ax=axlist[nv],bins=bins)
            axlist[nv].plot(times, times*0, c='k')
            axlist[nv].plot(times, times*0+1, c='k')
            axlist[nv].plot(times, times*0-1, c='k')
            axbonk(axlist[nv],xlabel=r'$t/t_{ff}$', ylabel=label)
            if symlog:
                axlist[nv].set_yscale('symlog',linthresh=1)

        outname='/home/dccollins/PigPen/%s_vr_vt_c%04d.png'%(this_looper.sim_name,core_id)
        fig.savefig(outname)
        print(outname)

sims=['u501', 'u502','u503']
for sim in sims:
    vr_vt(TL.loops[sim],symlog=True, core_list=[149,323])
    break

