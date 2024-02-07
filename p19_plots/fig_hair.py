
from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)


def buddy_hair(this_looper,core_list=None, suffix='', color_dict={}, what_to_do='centroid', external_ax=None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    if external_ax is None:
        fig,ax=plt.subplots(1,1)
    else:
        ax = external_ax
    x_ext=extents()
    y_ext=extents()
    z_ext=extents()
    mini_scrubbers={}
    times = thtr.times
    times.shape = times.size,1
    dont_axbonk=False
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

            if what_to_do == 'hair':
                ax.scatter( p[x][0,:].flatten(),p[y][0,:].flatten(),c=[color_base]*nparticles,s=0.1)
                ax.scatter( p[x][-1,:].flatten(),p[y][-1,:].flatten(),c='r',s=0.1)
                ax.plot( p[x], p[y], c=color_alpha, linewidth=0.1)
            elif what_to_do ==  'centroid':
                ax.plot( pmean[x], pmean[y], c=color_base)
            else:
                ax.plot( times, p[0], c=color_alpha)
                dont_axbonk=True


            h_ext = [x_ext,y_ext,z_ext][x]
            v_ext = [x_ext,y_ext,z_ext][y]

            xmin,xmax = h_ext.minmax+nar([-0.01,0.01])
            ymin,ymax = v_ext.minmax+nar([-0.01,0.01])
            
            text_x = pmean[x][0]
            text_y = pmean[y][0]
            #ax.text( text_x, text_y, r'$%s$'%core_id, color=color_dict[core_id])
            if 1:

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
                axbonk(ax,xlabel='x/L', ylabel='y/L', xlim=[xmin,xmax],ylim=[ymin,ymax])

        outname='plots_to_sort/%s_buddies_%s_%s.pdf'%(this_looper.sim_name,'xyz'[LOS],suffix)
        if external_ax is None:
            fig.savefig(outname)
            print(outname)
            plt.close(fig)




sims=['u501','u502','u503']

sims=['u502']
import find_other_cores
for sim in sims:
        subsets = [[214], [74], [112,113], [369]]
        dothis = [0,0,1,2]
        #subsets = [[369]]
        #dothis = [2]
        fig,axes=plt.subplots(2,2)
        this_looper=TL.loops[sim]

        for ns, core_list in enumerate(subsets):
            ax=axes.flatten()[ns]
            these_cores=core_list
            what_to_plot='hair'
            if dothis[ns] in [0]:
                color_dict={core_list[0]:[0.5]*4}
                suffix='c%04d'%core_list[0]
            elif dothis[ns] in [1]:
                b = [0, 90/256, 181/256, 0.5]
                r = [220/256, 50/256, 32/256, 0.5]

                color_dict = dict(zip(core_list,[r,b]))
                suffix='c%04d_c%04d'%tuple(core_list)
            else:
                what_to_plot='centroid'
                #these_cores=TL.loops[sim].core_by_mode['S2']
                #thtr=this_looper.tr
                #ms = trackage.mini_scrubber(thtr,369, do_velocity=False)
                #these_cores, shift = find_other_cores.get_other_cores( this_looper, [369], {369:ms})
                these_cores=[367,368,369,370,371]
                color_dict = colors.make_core_cmap(these_cores, cmap = 'autumn', seed = -1)
                #print(these_cores)
                #suffix='S2'
            buddy_hair(TL.loops[sim], core_list=these_cores,suffix=suffix,color_dict=color_dict, external_ax=ax, what_to_do=what_to_plot)
        fig.tight_layout()
        fig.savefig('plots_to_sort/F1_hair.pdf', bbox_inches='tight')

