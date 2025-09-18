
from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL
import movie_frames 

def mass_density(this_looper,core_list=None, do_plots=True, mass=None, dof=None, volume=None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    mask = movie_frames.quantized_mask(this_looper).flatten()
    times=thtr.times[mask]+0 #the zero makes a copy
    all_times=thtr.times
    times.shape=times.size,1
    times=times/colors.tff
    rho_all = thtr.track_dict['density']
    rho_min=rho_all.min()
    rho_max=rho_all.max()

    for nc,core_id in enumerate(core_list):
        print('V %s %d'%(this_looper.sim_name,core_id))
            
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True)
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

        dv = ms.cell_volume[sl].transpose()[mask,:]
        vv = dv.sum(axis=1)
        if do_plots:

            fig,axes=plt.subplots(2,1)
            ax=axes[0];ax1=axes[1]
            ax2=ax.twinx()

            ax.plot(times , rho, c=c, linewidth=0.1)
            axbonk(ax,xlabel=r'$t/t_{ff}$', ylabel=r'$\rho$',yscale='log', ylim=[rho_min,rho_max])

            ax2.plot( times, nar(mass[core_id])[mask],c='r')
            ax.plot( times, nar(dof[core_id])[mask],c='c')
            ax1.plot( times, nar(volume[core_id])[mask],c='g',label='Vtot')
            mean_volume = ms.cell_volume.mean(axis=0)
            ax1.plot( times, mean_volume[mask],label='<dv>')
            ax1.legend(loc=0)
            axbonk(ax1,yscale='log',xlabel='t/tff', ylabel='Volume')


            outname='plots_to_sort/%s_mass_volume_t_c%04d.png'%(this_looper.sim_name,core_id)
            fig.savefig(outname)



sims=['u501', 'u502','u503']
import three_loopers_u500 as TL
import mass_tools
if 'mt' not in dir():
    mt={}
for sim in sims:
    if sim not in mt:
        mt[sim]=mass_tools.mass_tool(TL.loops[sim])
        mt[sim].run()

sims=['u501', 'u502','u503']
for sim in sims:
    mass_density(TL.loops[sim], do_plots=True, mass=mt[sim].unique_mass, dof=mt[sim].dof, volume=mt[sim].volume)

