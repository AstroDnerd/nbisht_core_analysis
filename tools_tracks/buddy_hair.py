
from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)


def buddy_hair(this_looper,core_list=None, suffix='', color_dict={}, what_to_plot='centroids', shifter={}):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    fig,axlist=plt.subplots(1,3, figsize=(12,4))
    x_ext=extents()
    y_ext=extents()
    z_ext=extents()
    x_mean_ext=extents()
    y_mean_ext=extents()
    z_mean_ext=extents()
    mini_scrubbers={}
    times = (thtr.times+0)/colors.tff
    times.shape = times.size,1
    ds = this_looper.load(0)
    for core_id in core_list:
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
        ms.particle_pos(core_id)

        mini_scrubbers[core_id]=ms
        shift=np.zeros(3)
        if core_id in shifter:
            shift=shifter[core_id]

        x_ext( ms.particle_x+shift[0])
        y_ext( ms.particle_y+shift[1])
        z_ext( ms.particle_z+shift[2])
        x_mean_ext( ms.mean_x+ shift[0])
        y_mean_ext( ms.mean_y+ shift[1])
        z_mean_ext( ms.mean_z+ shift[2])
    for LOS in [0,1,2]:
        x = [1,0,1][LOS] # Using [1,0,1] and [2,2,0] 
        y = [2,2,0][LOS] # unfolds nicely.
        ax=axlist[LOS]
        #ax.set_aspect('equal')

        
        ycoord=[]
        for ncore,core_id in enumerate(core_list):
            ms=mini_scrubbers[core_id]
            

            if ms.nparticles < 1000:
                sl=slice(None)
                #c=[0.5]*4
                alpha=0.1
            else:
                sl = slice(None,None,10)
                #c=[0,0,0,0.1]
                #c=[0.1]*4
                alpha=0.1
            color_base = color_dict.get(core_id, [0.5]*4)[:3]
            color_alpha = np.concatenate([color_base,[alpha]])
            p = np.stack([ms.particle_x[sl].transpose(),ms.particle_y[sl].transpose(),ms.particle_z[sl].transpose()])
            pmean = np.stack([ms.mean_x, ms.mean_y, ms.mean_z])
            if core_id in shifter:
                sh = shifter[core_id]+0
                sh.shape = sh.size, 1,1
                p = p + sh
                sh.shape = sh.size,1
                pmean += sh
            nparticles = p[0][0,:].size

            dont_axbonk=False
            if what_to_plot == 'hair':
                ax.scatter( p[x][0,:].flatten(),p[y][0,:].flatten(),c=[color_base]*nparticles,s=0.1)
                ax.scatter( p[x][-1,:].flatten(),p[y][-1,:].flatten(),c='r',s=0.1)
                ax.plot( p[x], p[y], c=color_alpha, linewidth=0.1)
            elif what_to_plot == 'centroid':
                ax.plot( pmean[x], pmean[y], c=color_base)
            elif what_to_plot == 'xyzt':
                ax.plot( times, p[LOS], c=color_alpha, linewidth=0.1)
                ycoord.append(p[LOS].mean(axis=1)[-1])

                dont_axbonk=True


            if what_to_plot == 'centroid':
                h_ext = [x_mean_ext,y_mean_ext,z_mean_ext][x]
                v_ext = [x_mean_ext,y_mean_ext,z_mean_ext][y]
            else:
                h_ext = [x_ext,y_ext,z_ext][x]
                v_ext = [x_ext,y_ext,z_ext][y]
            xmin,xmax = h_ext.minmax+nar([0.0,0.05])
            ymin,ymax = v_ext.minmax#+nar([-0.01,0.01])
            
            #Text right by core.
            #text_x = pmean[x][0]
            #text_y = pmean[y][0]
            #ax.text( text_x, text_y, r'$%s$'%core_id, color=color_dict[core_id])

            if 0:

                this_p = [q[-1,:] for q in p]
                this_ax=ax

                the_x_tmp = this_p[x] - this_p[x].min()+ 0.1
                the_y_tmp = this_p[y] - this_p[y].min()+ 0.1
                the_x = 10**np.log10(the_x_tmp).mean() + this_p[x].min() - 0.1
                the_y = 10**np.log10(the_y_tmp).mean() + this_p[y].min() - 0.1
                text_height=( ymax-ymin)*0.05

                text_ymax = ymax-text_height
                text_y = text_ymax - ncore*text_height
                text_x = xmax
                this_ax.text( text_x, text_y, r'$%s$'%core_id, color=color_dict[core_id])
                this_ax.plot([the_x,text_x],[the_y,text_y],c='k',linewidth=0.1)
                particle_h = pmean[x][-1]
                particle_v = pmean[y][-1]
                text_height=0.02

        vert_order = np.argsort( ycoord)[::-1]
        for ncore,core_id in enumerate(nar(core_list)[vert_order]):
            ms=mini_scrubbers[core_id]
            ax.set(xlim=[times.min(),1.2*times.max()])

            if 1:
                #Use the Axis coordinates for the text.
                Hcoord=ds.coordinates.x_axis[LOS]
                Vcoord=ds.coordinates.y_axis[LOS]
                text_height=0.04
                start_y = 1-text_height
                start_x = 1-4*text_height
                text_y = start_y - ncore*text_height
                text_x = start_x
                points_axis=(text_x,text_y)
                axis_to_data = ax.transAxes + ax.transData.inverted()
                points_data = axis_to_data.transform(points_axis)
                ax.text( text_x,text_y,r'$c%04d$'%core_id, transform=ax.transAxes, color=color_dict.get(core_id,'k'))
            if 0:
                frame_ind=-1    
                if what_to_plot in ['hair','centroid']:
                    P = np.stack([ms.mean_x[frame_ind], ms.mean_y[frame_ind], ms.mean_z[frame_ind]])
                    #P += shifts[nc]
                    #ax.scatter( P[Hcoord], P[Vcoord])
                    ax.plot( [P[Hcoord], points_data[0]], [P[Vcoord], points_data[1]],linewidth=0.1, color=color_dict.get(core_id,'k'))
                else:
                    P = np.stack([ms.mean_x[frame_ind], ms.mean_y[frame_ind], ms.mean_z[frame_ind]])
                    #P += shifts[nc]
                    #ax.scatter( P[Hcoord], P[Vcoord])
                    color = color_dict.get(core_id,'k')
                    print(times.flatten()[-1], points_data, points_axis)
                    ax.plot( [times.flatten()[-1], points_data[0]], [P[LOS], points_data[1]], color=color)

        if not dont_axbonk:
            axbonk(ax,xlabel='xyz [code length]'[x], ylabel='xyz [code length]'[y], xlim=[xmin,xmax],ylim=[ymin,ymax])
        else:
            ax.set( xlabel=r'$t/t_{ff}$', ylabel = 'xyz'[LOS])

    outname='plots_to_sort/%s_buddies_%s_%s_%s.png'%(this_looper.sim_name, what_to_plot,'xyz',suffix)
    fig.tight_layout()
    fig.savefig(outname)
    print(outname)
    plt.close(fig)


