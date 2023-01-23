
from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL
import movie_frames 

def simple_g_dot_r(this_looper,core_list=None):

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
        ms.compute_ge(core_id)

        if ms.nparticles < 1000:
            sl=slice(None)
            c=[0.5]*4
        else:
            sl = slice(None,None,10)
            #c=[0,0,0,0.1]
            c=[0.1]*4

        g_dot_r = ms.rx_rel*ms.gx+ms.ry_rel*ms.gy+ms.rz_rel*ms.gz
        g2 = ms.gx*ms.gx+ms.gy*ms.gy+ms.gz*ms.gz
        g2 = g2[sl].transpose()[mask,:]
        r2 = ms.r2[sl].transpose()[mask,:]
        g_dot_r = g_dot_r[sl].transpose()[mask,:]
        #rho = ms.density[sl].transpose()
        #rho = rho[mask,:]

        ax.scatter(r2,g2)
        ax.set(xscale='log',yscale='log')

        if 0:
            cosgr= g_dot_r/np.sqrt(g2*r2)
            print(cosgr)
            pdb.set_trace()
            ax.plot(times ,cosgr, c=c, linewidth=0.1)
            #pdb.set_trace()
            axbonk(ax,xlabel=r'$t/t_{ff}$', ylabel=r'$g\cdot v$')#, ylim=[rho_min,rho_max])
            #ax.set_yscale('symlog',linthresh=10)
            import heat_map
            heat_map.heat_map(g_dot_r/np.sqrt(g2*r2),times.flatten(),ax)

        outname='plots_to_sort/%s_g_dot_r_t_c%04d.png'%(this_looper.sim_name,core_id)
        fig.savefig(outname)
        print(outname)




sims=['u501', 'u502','u503']
sims=['u502']
for sim in sims:
    core_list=[74]

    simple_g_dot_r(TL.loops[sim], core_list=core_list)

