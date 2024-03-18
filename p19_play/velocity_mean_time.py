
from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL
import movie_frames 

def vel_per_core(this_looper,core_list=None, tsing=None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    mask = movie_frames.quantized_mask(this_looper).flatten()
    times=thtr.times[mask]+0 #the zero makes a copy
    #times.shape=times.size,1
    times=times/colors.tff
    rho_all = thtr.track_dict['density']
    rho_min=rho_all.min()
    rho_max=rho_all.max()
    fig,ax=plt.subplots(3,1)
    ax0=ax[0];ax1=ax[1];ax2=ax[2]
    ax0.axhline(-1,c=[0.5]*4)
    ax1.axhline(1,c=[0.5]*4)
    ax2.axhline(1,c=[0.5]*4)
    for core_id in core_list:
        print('vel on core ', this_looper.sim_name,core_id)

            
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True)
        ms.particle_pos(core_id)

        my_times = times
        if tsing is not None:
            my_times = times / tsing.tsing_core[core_id]

        if True:
            sl=slice(None)
            c=[0.5]*4
        else:
            sl = slice(None,None,10)
            #c=[0,0,0,0.1]
            c=[0.1]*4

        rho = ms.density[sl].transpose()
        rho = rho[mask,:]

        v22 = ms.rel_vmag[sl].transpose()[mask,:]
        vr = ms.vr_rel[sl].transpose()[mask,:]
        vt = (ms.vt2_rel[sl].transpose()[mask,:])**0.5
        vrm=vr.mean(axis=1)
        v2 = v22.mean(axis=1)
        vtm=vt.mean(axis=1)

        ax0.plot( my_times, vrm, color=c)
        ax1.plot( my_times, vtm, color=c)
        ax2.plot( my_times, v2, color=c)




    outname='plots_to_sort/%s_mean_v_t.png'%(this_looper.sim_name)
    ax0.axhline(-1,c=c)
    fig.savefig(outname)
    print(outname)


import tsing
sim_list=['u501', 'u502','u503']
if 'tsing_tool' not in dir():
    tsing_tool={}
    for ns,sim in enumerate(sim_list):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()

sims=['u501', 'u502','u503']
usims=['u502']
for sim in sims:
    core_list = TL.loops[sim].core_by_mode['Alone']
    #core_list=[114]
    #core_list=list(core_list[:5])#+[114]   
    vel_per_core(TL.loops[sim], core_list=core_list, tsing=tsing_tool[sim])

