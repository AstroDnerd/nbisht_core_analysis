from starter2 import *
from collections import defaultdict
from scipy.optimize import curve_fit
import heat_map
reload(heat_map)
import movie_hair
reload(movie_hair)
import colors
import time_line
reload(time_line)
import movie_frames




if 1:
    import three_loopers_u500 as TL5
    this_looper=TL5.loops['u503']

if 1:
    core_list = [223]
    core_list = [76]
    full_core_list = np.unique( this_looper.tr.core_ids)
    core_list = full_core_list

if 0:
    movie_mask = movie_frames.quantized_mask(this_looper)
    times = this_looper.tr.times[movie_mask]/colors.tff
    frames = this_looper.tr.frames[movie_mask]

if 0:
    #hair movie
    #fig, ax = plt.subplots(1,1,figsize=(12,8))
    fig, ax = plt.subplots(1,1,figsize=(8,8))
    thing = movie_hair.flow(this_looper)
    thing.run(core_list=core_list, frames='reg', external_ax=ax, external_fig=fig)

if 1:
    #velocity heat map for each core.  
    reload(heat_map)
    #fig,ax=plt.subplots(1,1, figsize=(8,8))
    for field in ['velocity_x','velocity_y','velocity_z']:
        heat_map.heat_for_quantity( this_looper, field=field,core_list=core_list)#,external_ax=ax, bins=None)
        #axbonk(ax, xlabel=r'$t/t_{ff}$', ylabel=field, yscale='linear')
        #fig.savefig('plots_to_sort/heat_%s_%s_c%04d.png'%(field,this_looper.sim_name,core_list[0]))

    #name_template = 'plots_to_sort/heat_density_timeline_'+this_looper.sim_name+'_c%04d'%core_list[0]+'_n%04d'
    #time_line.timelines(fig,ax,times, frames=frames,name_template=name_template)

if 0:
    #Density heat map for each core.  
    reload(heat_map)
    fig,ax=plt.subplots(1,1, figsize=(8,8))
    heat_map.heat_for_quantity( this_looper, field='density',core_list=core_list,external_ax=ax, bins='plog')
    axbonk(ax, xlabel=r'$t/t_{ff}$', ylabel='density', yscale='log')
    #fig.savefig('plots_to_sort/heat_%s_%s_c%04d.png'%('density',this_looper.sim_name,core_list[0]))

    name_template = 'plots_to_sort/heat_density_timeline_'+this_looper.sim_name+'_c%04d'%core_list[0]+'_n%04d'
    time_line.timelines(fig,ax,times, frames=frames,name_template=name_template)
    
if 0:
    #potential  heat map for each core.  
    core_list=[76]
    full_core_list = np.unique( this_looper.tr.core_ids)
    core_list=full_core_list
    times = this_looper.tr.times/colors.tff

    for nc, core_id in enumerate(core_list):
        density = this_looper.tr.c([core_id],'PotentialField')
        cell_volume = this_looper.tr.c([core_id],'cell_volume')
        bins = np.linspace( density.min(), density.max(), 64)
        fig,ax=plt.subplots(1,1)
        heat_map.heat_map( density, times, bins=bins,ax=ax)

        mean_density=(density*cell_volume).sum(axis=0)/cell_volume.sum(axis=0)
        ax.plot( times, mean_density)
        outname='plots_to_sort/potential_time_heat_%s_c%04d.png'%(this_looper.sim_name,core_id)
        fig.savefig(outname)
        print(outname)
