
from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL
import movie_frames 

def rho_ke_ge(this_looper,core_list=None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    frames = movie_frames.quantized_mask(this_looper).flatten()
    times=thtr.times[frames]+0 #the zero makes a copy
    times.shape=times.size,1
    times=times/colors.tff
    rho_all = thtr.track_dict['density']
    rho_min=rho_all.min()
    rho_max=rho_all.max()
    for core_id in core_list:
        if 0:
            fig,ax=plt.subplots(3,1, figsize=(8,12))
            fig.subplots_adjust(hspace=0, wspace=0)
            ax0=ax[0]; ax1=ax[1]; ax2=ax[2]
        fig,ax=plt.subplots(1,1)

            
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True)
        ms.particle_pos(core_id)

        if ms.nparticles < 1000:
            sl=slice(None)
            c=[0.5]*4
        else:
            sl = slice(None,None,10)
            #c=[0,0,0,0.1]
            c=[0.1]*4

        phi = thtr.c([core_id],'PotentialField')[sl].transpose()
        phi = phi[frames,:]
        rho = ms.density[sl].transpose()
        rho = rho[frames,:]
        vx  = ms.raw_vx[sl].transpose()
        vx = vx[frames,:]
        vy  = ms.raw_vy[sl].transpose()
        vy = vy[frames,:]
        vz  = ms.raw_vz[sl].transpose()
        vz = vz[frames,:]


        v2  = vx**2 + vy**2 + vz**2
        ke  = 0.5*v2
        ge = 0.5*phi

        ax0=ax;ax1=ax;ax2=ax
        ax1=ax0.twinx()
        ax2=ax1
        c0 = [1.0,0.5,0.5,0.5]
        ax0.plot(times , rho, c=c0, linewidth=0.1)
        axbonk(ax0,xlabel=r'$t/t_{ff}$', ylabel=r'$\rho$',yscale='log', ylim=[rho_min,rho_max])

        c1 = [0.5,1.0,0.5,0.5]
        ax1.plot( times, ke, c=c1, linewidth=0.1)
        axbonk(ax1,xlabel=r'$t/t_{ff}$', ylabel=r'$KE$',yscale='log')

        c2 = [0.5,0.5,1.0,0.5]
        ax2.plot( times, np.abs(ge), c=c2, linewidth=0.1)
        axbonk(ax2,xlabel=r'$t/t_{ff}$', ylabel=r'$GE$',yscale='log')

        outname='plots_to_sort/%s_rho_ke_ge_t_c%04d.png'%(this_looper.sim_name,core_id)
        fig.savefig(outname)
        print(outname)




sims=['u501', 'u502','u503']
sims = ['u503']
for sim in sims:
    rho_ke_ge(TL.loops[sim], core_list=[189])

