
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
    G = colors.G
    gx = thtr.track_dict['grav_x']
    gy = thtr.track_dict['grav_y']
    gz = thtr.track_dict['grav_z']
    GE2 = -1/(8*np.pi)*(gx*gx+gy*gy+gz*gz)
    ge_min=GE2.min()
    ge_max=GE2.max()
    for core_id in core_list:
        fig,ax=plt.subplots(1,1)

            
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
        #ms.particle_pos(core_id)

        if ms.nparticles < 1000:
            sl=slice(None)
            c=[0.5]*4
        else:
            sl = slice(None,None,10)
            #c=[0,0,0,0.1]
            c=[0.1]*4

        rho = ms.density[sl].transpose()
        rho = rho[mask,:]
        gx = thtr.c([core_id],'grav_x')[sl].transpose()[mask,:]
        gy = thtr.c([core_id],'grav_y')[sl].transpose()[mask,:]
        gz = thtr.c([core_id],'grav_z')[sl].transpose()[mask,:]
        GE2 = 1/(8*np.pi*G)*(gx*gx+gy*gy+gz*gz)

        ax.plot(times , GE2, c=c, linewidth=0.1)
        axbonk(ax,xlabel=r'$t/t_{ff}$', ylabel=r'$(\nabla \phi)^2/8 pi G$',yscale='log', ylim=[ge_min,ge_max])
        ax2=ax.twinx()
        c=[1.0,0.1,0.1,0.1]
        ax2.plot(times , rho, c=c, linewidth=0.1)
        axbonk(ax2,xlabel=r'$t/t_{ff}$', ylabel=r'$\rho$',yscale='log')

        outname='plots_to_sort/%s_GE_t_c%04d.png'%(this_looper.sim_name,core_id)
        fig.savefig(outname)
        print(outname)




sims=['u501', 'u502','u503']
for sim in sims:
    core_list = np.unique(TL.loops[sim].tr.core_ids)
    core_list=[323]
    simple_rho(TL.loops[sim],core_list=core_list)
    break


