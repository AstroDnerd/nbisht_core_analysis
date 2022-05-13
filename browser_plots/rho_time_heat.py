
from starter2 import *
from collections import defaultdict
import scipy
import colors

import heat_map
reload(heat_map)

import three_loopers_u500 as TL
import movie_frames 

def simple_rho(this_looper,core_list=None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    mask = movie_frames.quantized_mask(this_looper).flatten()
    times=thtr.times[mask]/colors.tff
    rho_all = thtr.track_dict['density']
    rho_min=rho_all.min()
    rho_max=rho_all.max()
    for core_id in core_list:
        fig,ax=plt.subplots(1,1)

        quan = thtr.c([core_id],'density')[:,mask]
        bins = np.geomspace(rho_min,rho_max,64)
        XX,YY,HH,VV,plot=heat_map.heat_map( quan, times, ax=ax, bins=bins,hist_norm=True, zlim=[5e-3,3e-1])
        axbonk(ax,xlabel=r'$t/t_{ff}$', ylabel=r'$\rho/\rho_0$',yscale='log')
        norm = mpl.colors.LogNorm( HH[HH>0].min(), HH.max())
        fig.colorbar(plot)


        outname='plots_to_sort/%s_rho_t_heat_c%04d.png'%(this_looper.sim_name,core_id)
        fig.savefig(outname)
        print(outname)




sims=['u501', 'u502','u503']
for sim in sims:
    simple_rho(TL.loops[sim])#, core_list=[323])
#    break   

