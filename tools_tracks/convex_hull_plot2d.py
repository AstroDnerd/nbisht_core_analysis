

from starter2 import *
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
plt.close('all')
import time
from scipy.spatial import Delaunay
import figure_sublabel
import colors
import convex_hull_tools as CHT
def plot_2d(htool,core_list=None,accumulate=False,frames=[0],all_plots=False, label_cores=[],prefix="",
            color_dict=None,axis_to_plot=[0], plot_square=True, external_axis=None,
                plotstyle='xyz',
                add_jitter=True, center_image=True
           ):


    thtr = htool.this_looper.tr
    all_cores = np.unique(thtr.core_ids)
    rm = rainbow_map(len(all_cores))
    if frames is None:
        frames = thtr.frames
    if core_list is None:
        core_list = all_cores

    if external_axis is None:
        if -1 in axis_to_plot:
            fig_many, ax = plt.subplots(2,2,figsize=(8,8))
            ax=ax.flatten()
        else:
            fig_many, ax = plt.subplots(1,1)#,figsize=(4,4))
            ax=[ax] #because we need to refer to the first one.
            ax[0].set_aspect('equal')
    else:
        ax=external_axis
    if -1 in axis_to_plot:
        axistag='xyz'
    else:
        axistag='xyz'[axis_to_plot[0]]


    x_min, x_max, y_min, y_max = [1,0,1,0]
    x_ext = extents()
    y_ext = extents()
    z_ext = extents()
    if center_image:
        plot_x_ext = extents()
        plot_y_ext = extents()
    if color_dict is  None:
        color_dict = colors.make_core_cmap( core_list)

    ms_dict={}
    for ncore,core_id in enumerate(core_list):
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
        ms.particle_pos(core_id)
        ms_dict[core_id]=ms

        x_ext( ms_dict[core_id].particle_x)
        y_ext( ms_dict[core_id].particle_y)
        z_ext( ms_dict[core_id].particle_z)

    for it,frame in enumerate(frames):#asort):
        nt= np.where( nar(thtr.frames) == frame)[0][0]
        for aaa in ax:
            aaa.clear()
        for ncore,core_id in enumerate(core_list):
            ms = ms_dict[core_id]
            if ms.r.shape[0] <= 4:
                continue
            print('plot core %s %d'%(htool.this_looper.out_prefix,core_id))
            delta=0.1

            mask = slice(None)

            if not accumulate:
                if -1 in axis_to_plot:
                    for aaa in ax:
                        aaa.clear()
                        aaa.set_aspect('equal')
                        aaa.plot([0,1,1,0,0],[0,0,1,1,0])
                else:
                    ax.clear()
                    ax.set_aspect('equal')
                    ax.plot([0,1,1,0,0],[0,0,1,1,0])
            #this_x,this_y,this_z=ms.this_x[mask,nt],ms.this_y[mask,nt], ms.this_z[mask,nt]
            this_x,this_y,this_z=ms.particle_x[mask,nt],ms.particle_y[mask,nt], ms.particle_z[mask,nt]

            all_x,all_y,all_z=ms.particle_x,ms.particle_y, ms.particle_z
            all_p = [all_x,all_y,all_z]

            do_hull = True

            if np.unique(this_x).size < 4  or\
               np.unique(this_y).size < 4  or\
               np.unique(this_z).size < 4 :
                print("Not enough degrees of freedom")
                do_hull = False

            this_p = [this_x,this_y,this_z]


            if -1 in axis_to_plot:
                axis_to_actually_plot = [0,1,2]
            else:
                axis_to_actually_plot=axis_to_plot
            for LOS in axis_to_actually_plot:
                if plotstyle == 'xyz':
                    x = [1,0,1][LOS] # Using [1,0,1] and [2,2,0] 
                    y = [2,2,0][LOS] # unfolds nicely.
                if plotstyle == "Figure8":
                    x = [0,1,0][LOS] # Using [0,1,0] and [2,2,1] 
                    y = [2,2,1][LOS] # puts Y on the vert in panel 3
                xlab=r'$%s \rm(code\ length)$'%'xyz'[x]
                ylab=r'$%s \rm(code\ length)$'%'xyz'[y]


                if -1 in axis_to_plot or len(axis_to_plot)>1:
                    this_ax = ax[LOS]
                else:
                    this_ax = ax[0]
                n_particles = len(this_p[0])

                #
                # the plot
                #
                if add_jitter:
                    dx =(np.random.random( n_particles)-0.5)*1./128
                    dy =(np.random.random( n_particles)-0.5)*1./128
                else:
                    dx=0
                    dy=0
                if plotstyle=='Figure8':
                    size=3
                    marker=None
                else:
                    size=0.5
                    marker=None
                if center_image:
                    plot_x_ext( this_p[x]+dx)
                    plot_y_ext( this_p[y]+dy)
                this_ax.scatter(this_p[x]+dx, this_p[y]+dy,s=size, c=[color_dict[core_id]]*n_particles,
                                marker=marker, edgecolor='None')
                    

                #streaks.  Use with caution.
                #this_ax.plot(all_p[x].transpose(), all_p[y].transpose(), c=color_dict[core_id], linewidth=.1)

                if plot_square:
                    x_min = min([x_min,this_p[x].min(), -delta])
                    x_max = max([x_max,this_p[x].max(), 1+delta])
                    y_min = min([y_min,this_p[y].min(), -delta])
                    y_max = max([y_max,this_p[y].max(), 1+delta])

                    this_ax.plot([0,1,1,0,0], [0,0,1,1,0], c=[0.5]*3)
                else:
                    x_min, x_max = [x_ext,y_ext,z_ext][x].minmax
                    y_min, y_max = [x_ext,y_ext,z_ext][y].minmax


                if do_hull:
                    points_2d = np.array(list(zip(this_p[x],this_p[y])))
                    hull_2d = ConvexHull(points_2d)
                    vert_x = points_2d[hull_2d.vertices,0]
                    vert_y = points_2d[hull_2d.vertices,1]
                    vert_x = np.concatenate([vert_x,vert_x[0:1]])
                    vert_y = np.concatenate([vert_y,vert_y[0:1]])
                    this_ax.plot(vert_x, vert_y, 'k', linewidth=0.3)

                if core_id in label_cores or -1 in label_cores:
                    if plotstyle == 'Figure8':
                        if  LOS == 0:
                            #this_ax.text( this_p[x].max(), this_p[y].max(), r'$%s$'%core_id)
                            #the_x=this_p[x].mean()
                            #the_y=this_p[y].mean()
                            the_x_tmp = this_p[x] - this_p[x].min()+ 0.1
                            the_y_tmp = this_p[y] - this_p[y].min()+ 0.1
                            the_x = 10**np.log10(the_x_tmp).mean() + this_p[x].min() - 0.1
                            the_y = 10**np.log10(the_y_tmp).mean() + this_p[y].min() - 0.1
                            text_height=0.02
                            text_ymax = y_max-text_height
                            text_y = text_ymax - ncore*text_height
                            text_x = x_max - 2*text_height
                            this_ax.text( text_x, text_y, r'$%s$'%core_id, color=color_dict[core_id])
                            this_ax.plot([the_x,text_x],[the_y,text_y],c='k')
                    else:
                        text_x = this_p[x].mean()
                        text_y = this_p[y].max()
                        this_ax.text( text_x, text_y, r'$%s$'%core_id)



                def image_bumper(x_ext,y_ext, factor=1.1):
                    dx1 = max( [ np.abs(min([x_ext.minmax[0],0])), 
                                 np.abs(max([x_ext.minmax[1]-1,0]))])*factor
                    dy1 = max( [ np.abs(min([y_ext.minmax[0],0])), 
                                 np.abs(max([y_ext.minmax[1]-1,0]))])*factor
                    return max([dx1, dy1])
                if center_image:
                    ddd = image_bumper( plot_x_ext, plot_y_ext)
                    x_min, x_max=-ddd,1+ddd
                    y_min, y_max=-ddd,1+ddd
                #add sub labels.  Exits if theres nothing to do.
                figure_sublabel.add_sublabel(this_ax, htool.this_looper)
                axbonk(this_ax,xlabel=xlab,ylabel=ylab,xlim=[x_min,x_max],ylim=[y_min,y_max])
                #ax_many.set_title(title)
            cumltext=""
            if accumulate:
                cumltext="%04d"%ncore
            if all_plots and external_axis is None:
                outname = '%s/%s_hull_3d_t_%sc%04d_n%04d.png'%(dl.output_directory,htool.this_looper.out_prefix,cumltext,core_id,frame)
                fig_many.savefig(outname)
                print("Wrote "+outname)
        if 0:
            title_text=''
            if len(core_list) == 2:
                c1=core_list[0]; c2=core_list[1]
                n1=np.where(htool.cores_used==c1)[0][0]
                n2=np.where(htool.cores_used==c2)[0][0]
                o1 = htool.overlaps[c1][n2]
                o2 = htool.overlaps[c2][n1]
                ax[0].set_title('c %04d has %0.1e; c %04d has %0.1e'%(c1,o1, c2, o2))
            elif len(core_list) == 1:
                ax[0].set_title('c %04d is alone'%core_id)
            else:
                pass
        if accumulate and external_axis is None:
            outname = '%s/%s_hull_3d_t_%scXXXX_%s_n%04d.png'%(dl.output_directory,htool.this_looper.out_prefix,prefix,axistag,frame)
            fig_many.savefig(outname)
            print("Wrote "+outname)
