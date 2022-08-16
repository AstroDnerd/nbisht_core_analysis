
from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL
import movie_frames 

def simple_rho(this_looper,core_list=None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    mask = movie_frames.quantized_mask(this_looper).flatten()
    times=thtr.times[mask]+0 #the zero makes a copy
    times.shape=times.size,1
    times=times/colors.tff
    rho_all = thtr.track_dict['density']
    rho_min=rho_all.min()
    rho_max=rho_all.max()
    for core_id in core_list:
        fig,ax=plt.subplots(1,1)

            
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
        ms.particle_pos(core_id)

        if ms.nparticles < 1000:
            sl=slice(None)
            c=[0.5]*4
        else:
            sl = slice(None,None,10)
            #c=[0,0,0,0.1]
            c=[0.1]*4

        rho = ms.density[sl].transpose()
        rho = rho[mask,:]


        ax.plot(times , rho, c=c, linewidth=0.1)
        axbonk(ax,xlabel=r'$t/t_{ff}$', ylabel=r'$\rho$',yscale='log', ylim=[rho_min,rho_max])

        from scipy.ndimage import gaussian_filter

        UB = gaussian_filter(rho.max(axis=1),1)
        ax.plot( times, UB, c='k')
        ax2=plt.twinx()
        tf = times.flatten()
        dt = tf[1:]-tf[:-1]
        dU = UB[1:]-UB[:-1]
        tc = 0.5*(times[1:]+times[:-1])
        dU = dU[1:-1]
        dt = dt[1:-1]
        tc = tc[1:-1]
        dudt=dU/dt
        thresh = 1e5


        singularity = np.where( dudt >thresh)[0][0]

        if (dudt[singularity:]<=0).any():
            collapse_done = np.where( dudt[singularity:] < 0)[0][0]
            collapse_done += singularity
        else:
            collapse_done=-1
        ax2.plot(tc.flatten(), dudt, c='r')
        ax.scatter( tf[singularity], UB[singularity], marker='*', c='g',s=100)
        ax.scatter( tf[collapse_done], UB[collapse_done], marker='*',c='r')
        ax2.set_yscale('log')
        ax2.set_ylabel(r'$d\rho/dt$')
        ax2.set_ylim([rho_min,rho_max])
        ax2.axhline(1e5,c='r',linewidth=0.1)
        ax2.axvline(UB[singularity], c='k',linewidth=0.1)

        outname='plots_to_sort/%s_drho_dt_c%04d.png'%(this_looper.sim_name,core_id)
        fig.savefig(outname)
        print(outname)




sims=['u501', 'u502','u503']
#sims=['u502', 'u501']
for sim in sims:
    #core_list={'u501':[323], 'u502':[381], 'u503':[203]}[sim]
    core_list=None
    simple_rho(TL.loops[sim], core_list=core_list)

