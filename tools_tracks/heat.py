
from starter2 import *
from collections import defaultdict
import heat_map
from scipy.optimize import curve_fit

if 1:
    #mean density.  If you divide by the wrong volume you get a strange plot.
    
    import heat_map
    reload(heat_map)
    fig,ax=plt.subplots(1,1)
    this_looper=brthing.this_looper
    density = this_looper.tr.track_dict['density']
    minmin = density.min()
    maxmax = density.max()
    bins = np.geomspace(minmin,maxmax,64)
    heat_map.plot_heat(times = this_looper.tr.times, cores_used = brthing.cores_used, quan_dict=brthing.mean_rho ,ax=ax, cut_last=True,bins=bins)
    axbonk(ax,yscale='log', ylim=[minmin,maxmax])
    fig.savefig('plots_to_sort/mean_rho_per_core_%s.png'%brthing.this_looper.sim_name)


if 0:
    #Density heat map for each core.  Very cool.  Finish this.
    import colors
    core_list=[76]
    full_core_list = np.unique( this_looper.tr.core_ids)
    times = this_looper.tr.times/colors.tff

    for nc, core_id in enumerate(core_list):
        density = this_looper.tr.c([core_id],'density')
        cell_volume = this_looper.tr.c([core_id],'cell_volume')
        bins = np.geomspace( density.min(), density.max(), 64)
        fig,ax=plt.subplots(1,1)
        heat_map.heat_map( density, times, bins=bins,ax=ax)

        mean_density=(density*cell_volume).sum(axis=0)/cell_volume.sum(axis=0)
        ax.plot( times, mean_density)
        axbonk(ax,yscale='log')
        outname='plots_to_sort/density_time_heat_%s_c%04d.png'%(this_looper.sim_name,core_id)
        fig.savefig(outname)
        print(outname)
    plt.close('all')
    fig,ax=plt.subplots(1,1)
    bins = np.geomspace(mean_density_all.min(), mean_density_all.max(),64)
    heat_map.heat_map( mean_density_all,times, bins=bins,ax=ax)
    fig.savefig('plots_to_sort/mean_density_%s.png'%(this_looper.sim_name))

