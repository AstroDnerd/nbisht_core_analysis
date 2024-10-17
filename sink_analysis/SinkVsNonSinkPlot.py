import warnings
warnings.filterwarnings("ignore")
import matplotlib.ticker as tck
from matplotlib.ticker import FuncFormatter, MultipleLocator
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib import rc
rc('text', usetex=False) # Use LaTeX font
import h5py
import numpy as np
import random 
from starter2 import *
import track_loader as TL
import json
import pandas as pd

#which_cores = ['all', 'densest']
def SinkVsNonSinkPlotter(trackname, sink_trackname, which_cores = 'densest'):
    make_dir("plots_to_sort/SinkVsNonSink")
    this_track = track_info.tracks[trackname]
    this_simname = this_track.sim_directory
    with open(this_track.SinkClumpLink_fname, 'r') as fp:
        sink_core_dict = json.load(fp)
    TL.load_tracks([trackname])
    this_looper = TL.tracks[trackname]
    thtr=this_looper.tr
    main_core_list = nar(this_looper.core_ids)
    main_core_list.sort()
    tracklist = this_looper.frame_list
    mini_scrubbers={}
    for core_id in main_core_list:
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
        ms.particle_pos(core_id)
        mini_scrubbers[core_id]=ms
    if which_cores == 'densest':
        #get cores with total density above the mean total density of all cores
        all_dens = []
        for key_core in main_core_list:
            all_dens.append(mini_scrubbers[key_core].density_tot[-1])
        mean_dens = sum(all_dens)/len(all_dens)
        print("mean_dens of cores")
        new_core_list = []
        for key_id in range(len(main_core_list)):
            if all_dens[key_id]>mean_dens:
                new_core_list.append(main_core_list[key_id])
        main_core_list = new_core_list
        print("Plotting densest cores")
    elif which_cores == 'all':
        print("Plotting all cores")
    else:
        print("Invalid which_cores")
        return
    print("Number of cores plotting: ",len(main_core_list))

    for this_frame in tracklist[-1:]:
        ds = yt.load("%s/DD%04d/data%04d"%(this_simname,this_frame,this_frame))
        all_data = ds.all_data()

        #load sink simulation
        sink_track = track_info.tracks[sink_trackname]
        sink_simname = sink_track.sim_directory
        sink_frame = sink_track.target_frame
        got_right_frame = False
        #make sure sink frame chosen is same or close to non sink chosen frame
        while got_right_frame==False:
            ds_sink = yt.load("%s/DD%04d/data%04d"%(sink_simname,sink_frame,sink_frame))
            if ds_sink.current_time-ds.current_time >1e-6:
                sink_frame-=1
            elif ds_sink.current_time-ds.current_time <-1e-6:
                sink_frame+=1
            else:
                got_right_frame = True
        all_data_sink = ds_sink.all_data()
        for core_id in main_core_list:
            ms = mini_scrubbers[core_id]
            print('Frame %04d, core %03d'%(this_frame,core_id))
            ms.particle_pos(core_id)

            pdb.set_trace()
            #setup plot
            plot = yt.ProjectionPlot(ds, 'z', ("enzo","Density"), fontsize = 10)
            plot.annotate_timestamp(corner="upper_left", redshift=True, draw_inset_box=True)
            plot.annotate_scale(corner="upper_right")
            suppress_annotation = False

            #setup plot
            fig = plt.figure()

            # See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
            grid = AxesGrid(
                fig,
                (0.085, 0.085, 0.83, 0.83),
                nrows_ncols=(3, 2),
                axes_pad=0.05,
                label_mode="L",
                share_all=True,
                cbar_location="right",
                cbar_mode="edge",
                cbar_size="3%",
                cbar_pad="0%",
                aspect=False,
            )

            weight_arr = ["cell_volume", "cell_mass", "velocity_magnitude"]
            weight_lim = [[1e-10,1e-6],[1e-11,1e-1],[1e-1,1e+2]]
            i = 0
            j=0
            z_data_arr = []
            for Brho_weight in weight_arr:   
                p1 = yt.PhasePlot(ds,("gas","density"),("gas","magnetic_field_strength"), ("gas",Brho_weight), fontsize = 10)
                p2 = yt.PhasePlot(ds_sink,("gas","density"),("gas","magnetic_field_strength"), ("gas",Brho_weight), fontsize = 10)
                for p in [p1,p2]:
                    # Ensure the axes and colorbar limits match for all plots
                    p.set_xlim(1.0e-4, 1.0e+9)
                    p.set_ylim(1.0e-1, 1.0e+5)
                    plot = p.plots[("gas",Brho_weight)]
                    plot.figure = fig
                    plot.axes = grid[i].axes
                    if i%2==0:
                        p.set_zlim(("gas", Brho_weight), weight_lim[j][0], weight_lim[j][1])
                        plot.cax = grid.cbar_axes[j]

                    # Actually redraws the plot.
                    p.render()

                    # Modify the axes properties **after** p.render() so that they
                    # are not overwritten.
                    plot.axes.xaxis.set_minor_locator(plt.LogLocator(base=10.0, subs=[2.0, 5.0, 8.0]))
                    #make isocontours
                    xy =  np.array(plot.image.get_coordinates())
                    x_quad = xy[0,:,0]
                    y_quad = xy[:,0,1]
                    x_cen = (x_quad[:-1] + x_quad[1:]) / 2
                    y_cen = (y_quad[:-1] + y_quad[1:]) / 2

                    z_data = np.nan_to_num(np.array(plot.image.get_array()))
                    z_data_arr.append(z_data)
                    plot.axes.contour(x_cen,y_cen,z_data)
                    if p==p2:
                        norm = mpl.colors.Normalize(vmin=weight_lim[j][0], vmax=weight_lim[j][1])
                        for sink_i in range(len(sink_particle_mass[-1])):
                            sink_pos = [float(sink_particle_xpos[-1][sink_i]), float(sink_particle_ypos[-1][sink_i]),float(sink_particle_zpos[-1][sink_i])]
                            point_obj =  ds_sink.point(sink_pos)
                            sink_cell_volume = point_obj['gas', 'cell_volume']
                            sink_cell_density = sink_particle_mass[-1][sink_i]/float(sink_cell_volume)
                            sink_cell_magfield = float(point_obj['gas', 'magnetic_field_strength'])
                            color_variable_val = float(point_obj['gas', Brho_weight])
                            c_sink = plot.cb.cmap(norm(color_variable_val))
                            plot.axes.scatter(sink_cell_density,sink_cell_magfield, s = 5, c = c_sink)

                    i+=1
                j+=1



SinkVsNonSinkPlotter('nb101', 'nb102')