
from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL

def grader(this_looper,core_list=None, mode='One', mini_scrubbers=None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)
    all_cores=np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    times = thtr.times
    times.shape = times.size,1
    if mini_scrubbers is None:

        mini_scrubbers={}
        for core_id in all_cores:
            print('ms',core_id)
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            ms.particle_pos(core_id)

            mini_scrubbers[core_id]=ms

    other_cores={}
    for main_core in core_list:
        fig,axes=plt.subplots(2,2)
        fig.subplots_adjust(hspace=0,wspace=0)

        x_ext=extents()
        y_ext=extents()
        z_ext=extents()
        other_cores[main_core]=[]
        ms = mini_scrubbers[main_core]
        mmx = ms.mean_x.mean()
        mmy = ms.mean_y.mean()
        mmz = ms.mean_z.mean()
        thresh = 0.05
        left =  nar([mmx-thresh,mmy-thresh,mmz-thresh])
        right = nar([mmx+thresh,mmy+thresh,mmz+thresh])
        x_ext(nar([mmx-thresh,mmx+thresh]))
        y_ext(nar([mmy-thresh,mmy+thresh]))
        z_ext(nar([mmz-thresh,mmz+thresh]))
        x_ext(ms.mean_x)
        y_ext(ms.mean_y)
        z_ext(ms.mean_z)
        for ocore in all_cores:
            if ocore==main_core:
                continue
            oms = mini_scrubbers[ocore]
            if 0:
                ok = ((oms.particle_x > left[0])*\
                        (oms.particle_x < right[0])*\
                        (oms.particle_y > left[1])*\
                        (oms.particle_y < right[1])*\
                        (oms.particle_z > left[2])*\
                        (oms.particle_z < right[2])).any()
            ok = ( ((oms.mean_x-ms.mean_x)**2+(oms.mean_y-ms.mean_y)**2+(oms.mean_z-ms.mean_z)**2) < 0.01).any()
            if ok:
                other_cores[main_core].append(ocore)
                x_ext(oms.mean_x)
                y_ext(oms.mean_y)
                z_ext(oms.mean_z)

        color_dict = colors.make_core_cmap(other_cores[main_core], cmap = 'autumn', seed = -1)
        color_dict[main_core]=[0.0,1.0,0.0]
        cores_to_plot = other_cores[main_core]+[main_core]
        print(cores_to_plot)
        for LOS in [0,1,2]:
            ax=axes.flatten()[LOS]
            x = [1,0,1][LOS] # Using [1,0,1] and [2,2,0] 
            y = [2,2,0][LOS] # unfolds nicely.
            for ncore,core_id in enumerate(cores_to_plot):
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
                if 0:
                    ax.scatter( p[x][0,:].flatten(),p[y][0,:].flatten(),c=[color_base]*nparticles,s=0.1)
                    ax.scatter( p[x][-1,:].flatten(),p[y][-1,:].flatten(),c='r',s=0.1)
                    ax.plot( p[x], p[y], c=color_alpha, linewidth=0.1)
                elif 1:
                    ax.plot( pmean[x], pmean[y], c=color_base)
                else:
                    ax.plot( times, p[0], c=color_alpha)
                    dont_axbonk=True


                #xmin,xmax = left[x],right[x]
                #ymin,ymax = left[y],right[y]
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

        axes[0][0].set(xticks=[],xlabel="")
        axes[0][1].set(yticks=[],ylabel="",xticks=[],xlabel="")
        #axes[1][0].set(yticks=[],ylabel="",xticks=[],xlabel="")
        axes[1][1].set(yticks=[],ylabel="")
        outname='plots_to_sort/%s_grader_c%04d.png'%(this_looper.sim_name,main_core)
        fig.savefig(outname)
        print(outname)
        plt.close(fig)


if 'mini_scrubbers' not in dir():
    mini_scrubbers={}
    for sim in ['u501','u502','u503']:
        this_looper=TL.loops[sim]
        mini_scrubbers[sim]={}
        all_cores=np.unique(this_looper.tr.core_ids)
        thtr=this_looper.tr

        for core_id in all_cores:
            print('ms',core_id)
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            ms.particle_pos(core_id)

            mini_scrubbers[sim][core_id]=ms

#sims=['u501']#, 'u502','u503']
#sims=['u501','u502','u503']
sims=['u501']
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

    
S4 = [ 241, 235, 236, 237, 212, 228, 233, 246, 232, 234, 230, 222, 226]
S5 = [128, 146, 127, 144, 142, 126, 104, 103]
S6 = [93,94,95,96,97,98,99]

for sim in sims:
    mode_list=TL.loops[sim].unique_modes
    #core_list=[295]
    core_list=[8]
    core_list=None
    grader(TL.loops[sim], core_list=core_list, mode=mode, mini_scrubbers=mini_scrubbers[sim])
