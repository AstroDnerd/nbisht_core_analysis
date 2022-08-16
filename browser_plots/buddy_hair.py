
from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL

def buddy_hair(this_looper,core_list=None, mode='One', color_dict={}, what_to_plot='centroids'):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    fig,ax=plt.subplots(1,1)
    x_ext=extents()
    y_ext=extents()
    z_ext=extents()
    mini_scrubbers={}
    times = thtr.times
    times.shape = times.size,1
    for core_id in core_list:
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
        ms.particle_pos(core_id)

        mini_scrubbers[core_id]=ms
        x_ext( ms.particle_x)
        y_ext( ms.particle_y)
        z_ext( ms.particle_z)
    for LOS in [2]:
        x = [1,0,1][LOS] # Using [1,0,1] and [2,2,0] 
        y = [2,2,0][LOS] # unfolds nicely.
        for ncore,core_id in enumerate(core_list):
            ms=mini_scrubbers[core_id]
            

            if ms.nparticles < 1000:
                sl=slice(None)
                #c=[0.5]*4
                alpha=0.5
            else:
                sl = slice(None,None,10)
                #c=[0,0,0,0.1]
                #c=[0.1]*4
                alpha=0.1
            color_base = color_dict.get(core_id, [0.5]*4)[:3]
            color_alpha = np.concatenate([color_base,[alpha]])
            p = [ms.particle_x[sl].transpose(),ms.particle_y[sl].transpose(),ms.particle_z[sl].transpose()]
            pmean = [ms.mean_x, ms.mean_y, ms.mean_z]
            nparticles = p[0][0,:].size

            dont_axbonk=False
            if what_to_plot == 'hair':
                ax.scatter( p[x][0,:].flatten(),p[y][0,:].flatten(),c=[color_base]*nparticles,s=0.1)
                ax.scatter( p[x][-1,:].flatten(),p[y][-1,:].flatten(),c='r',s=0.1)
                ax.plot( p[x], p[y], c=color_alpha, linewidth=0.1)
            elif what_to_plot == 'centroids':
                ax.plot( pmean[x], pmean[y], c=color_base)
            else:
                ax.plot( times, p[0], c=color_alpha)
                dont_axbonk=True


            h_ext = [x_ext,y_ext,z_ext][x]
            v_ext = [x_ext,y_ext,z_ext][y]
            xmin,xmax = h_ext.minmax+nar([-0.1,0.1])
            ymin,ymax = v_ext.minmax+nar([-0.1,0.1])
            
            text_x = pmean[x][0]
            text_y = pmean[y][0]
            ax.text( text_x, text_y, r'$%s$'%core_id, color=color_dict[core_id])
            if 0:

                this_p = [q[-1,:] for q in p]
                this_ax=ax

                the_x_tmp = this_p[x] - this_p[x].min()+ 0.1
                the_y_tmp = this_p[y] - this_p[y].min()+ 0.1
                the_x = 10**np.log10(the_x_tmp).mean() + this_p[x].min() - 0.1
                the_y = 10**np.log10(the_y_tmp).mean() + this_p[y].min() - 0.1
                text_height=0.04
                text_ymax = ymax-text_height
                text_y = text_ymax - ncore*text_height
                text_x = xmax - 2*text_height
                this_ax.text( text_x, text_y, r'$%s$'%core_id, color=color_dict[core_id])
                this_ax.plot([the_x,text_x],[the_y,text_y],c='k',linewidth=0.1)
                particle_h = pmean[x][-1]
                particle_v = pmean[y][-1]
                text_height=0.02


        if not dont_axbonk:
            axbonk(ax,xlabel='xyz [code length]'[x], ylabel='xyz [code length]'[y], xlim=[xmin,xmax],ylim=[ymin,ymax])

        outname='plots_to_sort/%s_buddies__%s_%s_%s.pdf'%(this_looper.sim_name, what_to_plot,'xyz'[LOS],mode)
        fig.savefig(outname)
        print(outname)
        plt.close(fig)




#sims=['u501']#, 'u502','u503']
sims=['u501','u502','u503']
#sims=['u503']
if 0:
    #get them all to sanitize
    all_modes = []
    for sim in sims:
        all_modes += TL.loops[sim].unique_modes
    all_modes=np.unique(nar(sorted(all_modes)))

    print(all_modes)
#Skip mostly because they're not fleshed out.
#No cores are missed by this cut.
skip = [ 'Barely', 'Found', 'Merger', 'Odd',  'Shard', 'Tides']
#this one skips singletons.
skip += [ 'One']

for sim in sims:
    mode_list=TL.loops[sim].unique_modes
    #mode_list=['S4']
    for mode in mode_list:
        if mode in skip:
            continue
        core_list = set(TL.loops[sim].core_by_mode[mode])
        #core_list = list(core_list - set(S4)-set(S5) - set(S6))

        color_dict = colors.make_core_cmap(core_list, cmap = 'tab20c', seed = -1)
        #print(mode, core_list)
        buddy_hair(TL.loops[sim], core_list=core_list, mode=mode,color_dict=color_dict, what_to_plot='centroid')

