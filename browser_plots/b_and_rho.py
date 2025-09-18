
from starter2 import *
from collections import defaultdict
import scipy
import colors
reload(colors)

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
    B_all = thtr.track_dict['magnetic_field_strength']
    B_min=B_all.min()
    B_max=B_all.max()
    for core_id in core_list:
        fig,ax_square=plt.subplots(2,2)
        ax=ax_square.flatten()
        

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
        B=thtr.c([core_id],'magnetic_field_strength')[sl].transpose()[mask,:]#/colors.mean_field[this_looper.sim_name]

        ax[0].plot(times , rho, c=c, linewidth=0.1)
        axbonk(ax[0],xlabel=r'$t/t_{ff}$', ylabel=r'$\rho$',yscale='log', ylim=[rho_min,rho_max])

        ax[1].plot(times, B, c=c,linewidth=0.1)
        axbonk(ax[1], xlabel=r'$t/t_{ff}$',ylabel='B', yscale='log')

        ax[2].plot(times, np.log10(B)/np.log10(rho), c=c,linewidth=0.1)
        ax[2].set_yscale('symlog',linthresh=1)

        ax[3].scatter(rho[0,:], B[0,:])
        ax[3].scatter(rho[-1,:], B[-1,:])
        ax[3].plot( rho, B, c=c,linewidth=0.1)
        axbonk(ax[3],xscale='log',yscale='log',xlabel='rho',ylabel='B')

        outname='plots_to_sort/%s_b_and_rho_t_c%04d.png'%(this_looper.sim_name,core_id)
        fig.savefig(outname)
        print(outname)




sims=['u501', 'u502','u503']
for sim in sims:
    simple_rho(TL.loops[sim])#, core_list=[323])

