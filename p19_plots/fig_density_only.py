
from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)

import track_loader as TL

import movie_frames 

def simple_rho(this_looper,core_list=None, thicker=False, tsing=None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    if len(core_list) > 10:
        print('Cannot plot so many cores.')
        return
    thtr=this_looper.tr
    mask = movie_frames.quantized_mask(this_looper).flatten()
    times=thtr.times[mask]+0 #the zero makes a copy
    times.shape=times.size,1
    times=times/colors.tff
    rho_all = thtr.track_dict['density']
    rho_min=rho_all.min()*colors.density_units
    rho_max=rho_all.max()*colors.density_units
    fig,axes=plt.subplots(len(core_list),1, figsize=(4,4))
    fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    for nc,core_id in enumerate(core_list):
        print('eat ',this_looper.sim_name,core_id)
        ax=axes[nc]

            
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
        ms.particle_pos(core_id)

        if ms.nparticles < 1000:
            sl=slice(None,None,5)
            c=[0.5]*4
        else:
            sl = slice(None,None,15)
            #c=[0,0,0,0.1]
            c=[0.1]*4
        #sl = slice(None,None,5)
        c=[0.1]*4
        if thicker:
            c=[0.5]*3

        rho = ms.density[sl].transpose()
        rho = rho[mask,:]
        rho *= colors.density_units

        sim_number=int(this_looper.sim_name[-1])
        #ax.text(0,1e6,'sim%d core %d %s'%(sim_number, core_id, this_looper.mode_dict[core_id][-1]))
        ax.text(0,1e6,'sim%d core %d'%(sim_number, core_id))
        ax.plot(times , rho, c=c, linewidth=0.1)
        if tsing is not None:
            ax.axvline( tsing.tsing_core[core_id], c=[0.5]*3,linewidth=1)
            ax.axvline( tsing.tend_core[core_id], c=[0.5]*3,linewidth=1)
        axbonk(ax,xlabel=r'$t/t_{ff}$', ylabel=r'$n~[\rm{cm^{-3}]}$',yscale='log', ylim=[rho_min,rho_max])

    for na,aa in enumerate(axes):
        if na < len(axes)-1:
            aa.set(xticks=[])

    outname='plots_to_sort/%s_rho_t_several.pdf'%(this_looper.sim_name)
    print(outname)
    fig.savefig(outname, bbox_inches='tight')



import tsing
reload(tsing)
sim_list = ['u501','u502','u503']
TL.load_tracks(sim_list)
if 'tsing_tool' not in dir() or True:
    tsing_tool={}
    for ns,sim in enumerate(sim_list):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()


#sims=['u501', 'u502','u503']
sims=['u502']
for sim in sims:
    core_list = [214, 74]#, 112,113, 368]
    simple_rho(TL.loops[sim],core_list=core_list, thicker=True,tsing=tsing_tool[sim])

