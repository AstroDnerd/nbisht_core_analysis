
from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL

def simple_hair(this_looper,core_list=None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    for core_id in core_list:
        fig,ax=plt.subplots(1,1)

            
        for LOS in [0]:
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            ms.particle_pos(core_id)
            x = [1,0,1][LOS] # Using [1,0,1] and [2,2,0] 
            y = [2,2,0][LOS] # unfolds nicely.

            if ms.nparticles < 1000:
                sl=slice(None)
                c=[0.5]*4
            else:
                sl = slice(None,None,10)
                #c=[0,0,0,0.1]
                c=[0.1]*4
            p = [ms.particle_x[sl].transpose(),ms.particle_y[sl].transpose(),ms.particle_z[sl].transpose()]

            ax.scatter( p[x][0,:].flatten(),p[y][0,:].flatten(),c='k',s=0.1)
            ax.scatter( p[x][-1,:].flatten(),p[y][-1,:].flatten(),c='r',s=0.1)
            ax.plot( p[x], p[y], c=c, linewidth=0.1)
            axbonk(ax,xlabel='xyz [code length]'[x], ylabel='xyz [code length]'[y])

        fig.savefig('plots_to_sort/%s_hair_%s_c%04d.png'%(this_looper.sim_name,'xyz'[LOS],core_id))
        print('plots_to_sort/%s_hair_%s_c%04d.png'%(this_looper.sim_name,'xyz'[LOS],core_id))




sims=['u501', 'u502','u503']
for sim in sims:
    simple_hair(TL.loops[sim])

