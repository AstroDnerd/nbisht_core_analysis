import warnings
warnings.filterwarnings("ignore")
import matplotlib.ticker as tck
from matplotlib.ticker import FuncFormatter, MultipleLocator
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid
rc('text', usetex=False) # Use LaTeX font
import h5py
import numpy as np
import random 
from starter2 import *
import track_loader as TL
import json
import pandas as pd
import animate_plots

#which_cores = ['all', 'densest', 'alone]
#suppress_annotation for printing mass of sinks
#set_clim across all timeseries, basically does the loop twice, may take quite long lol, ignores clim_arr
def SinkVsNonSinkPlotter(trackname, sink_trackname, field_name = ("gas","density"), which_cores = 'densest', suppress_annotation = False, set_clim = True, clim_arr = [0,0]):
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
    elif which_cores == 'alone':
        new_core_list = []
        mode_file = h5py.File(this_track.mode_fname, 'r')
        all_cores = mode_file['modes'][()]
        for core_id in all_cores:
            if str(core_id[0])[2] == 'A' and int(str(core_id[0])[3:-1]) in main_core_list:
                new_core_list.append(int(str(core_id[0])[3:-1]))
        mode_file.close()
        main_core_list = new_core_list
        print("Plotting alone cores")
    else:
        print("Invalid which_cores")
        return
    print("Number of cores plotting: ",len(main_core_list))
    print("Cores: ",main_core_list)
    max_val = -1
    min_val = 1e64
    for this_frame_index in range(len(tracklist)):
        this_frame = tracklist[this_frame_index]
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
        if set_clim==True:
            max_ds = float(np.max(all_data[field_name]))
            min_ds = float(np.min(all_data[field_name]))
            max_ds_sink = float(np.max(all_data_sink[field_name]))
            min_ds_sink = float(np.min(all_data_sink[field_name]))
            max_val = max(max_ds/10,max_ds_sink/10,max_val)
            min_val = min(min_ds,min_ds_sink,min_val)
        solar_mass = 5900*1.91e33
        sink_particle_mass = []
        sink_particle_xpos = []
        sink_particle_ypos = []
        sink_particle_zpos = []
        number_of_sink_particles = []
        sink_particle_index = []
        time_arr = []

        type_4_indices = np.where(np.array(all_data_sink[ds_sink.fields.all.particle_type])==4)[0]
        sink_particle_mass.append(np.array(all_data_sink[ds_sink.fields.all.particle_mass])[type_4_indices].tolist())
        sink_particle_xpos.append(np.array(all_data_sink[ds_sink.fields.all.particle_position_x])[type_4_indices].tolist())
        sink_particle_ypos.append(np.array(all_data_sink[ds_sink.fields.all.particle_position_y])[type_4_indices].tolist())
        sink_particle_zpos.append(np.array(all_data_sink[ds_sink.fields.all.particle_position_z])[type_4_indices].tolist())
        number_of_sink_particles.append(len(type_4_indices))
        sink_particle_index.append(np.array(all_data_sink[ds_sink.fields.all.particle_index])[type_4_indices].tolist())
        time_arr.append(ds.current_time)

        for core_id in main_core_list:
            make_dir("plots_to_sort/SinkVsNonSink_%s_%s/Core_%03d"%(field_name[1],which_cores,core_id))
            ms = mini_scrubbers[core_id]
            print('Frame %04d, core %03d'%(this_frame,core_id))
            ms.particle_pos(core_id)
            core_center = ms.mean_center[:,this_frame_index]
            core_rms = ms.rms[this_frame_index]
            core_radius = core_rms+0.01
            ds_sphere = ds.sphere(core_center, core_radius)
            ds_sink_sphere = ds_sink.sphere(core_center, core_radius)
            #setup plot
            fig = plt.figure()
            # See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
            grid = AxesGrid(
                fig,
                (0.055, 0.055, 0.85, 0.85),
                nrows_ncols=(2, 3),
                axes_pad=0.8,
                label_mode="all",
                share_all=True,
                cbar_location="right",
                cbar_mode="edge",
                cbar_size="3%",
                cbar_pad="0%",
                aspect=False,
            )

            dim_arr = ['x','y','z']
        
            i=0
            z_data_arr = []
            for dim in dim_arr:   
                p1 = yt.ProjectionPlot(ds, dim, field_name, fontsize = 10, data_source =ds_sphere, center = core_center, width = 2*core_radius)
                p1.annotate_timestamp(corner="upper_left", draw_inset_box=True)
                p1.annotate_scale(corner="upper_right")
                if set_clim == False:
                    p1.set_zlim("all",clim_arr[0], clim_arr[1])
                
                p2 = yt.ProjectionPlot(ds_sink, dim, field_name, fontsize = 10, data_source =ds_sink_sphere, center = core_center, width = 2*core_radius)
                p2.annotate_timestamp(corner="upper_left", draw_inset_box=True)
                p2.annotate_scale(corner="upper_right")
                if set_clim == False:
                    p2.set_zlim("all",clim_arr[0], clim_arr[1])
                for sink_i in range(len(sink_particle_mass[-1])):
                    new_x = float(sink_particle_xpos[-1][sink_i]) - core_center[0]
                    new_y = float(sink_particle_ypos[-1][sink_i]) - core_center[1]
                    new_z = float(sink_particle_zpos[-1][sink_i]) - core_center[2]
                    if (new_x<core_radius and new_x>-core_radius) and (new_y<core_radius and new_y>-core_radius) and (new_z<core_radius and new_z>-core_radius):
                        p2.annotate_text((float(sink_particle_xpos[-1][sink_i]), float(sink_particle_ypos[-1][sink_i]), float(sink_particle_zpos[-1][sink_i])),
                                         '.', coord_system="data", text_args={"fontsize": 100, "color":'black'})
                        if suppress_annotation ==False:
                            offset = 0.001
                            p2.annotate_text((float(sink_particle_xpos[-1][sink_i])+offset, float(sink_particle_ypos[-1][sink_i])+offset, float(sink_particle_zpos[-1][sink_i])+offset),
                                             '{:0.1e}'.format(sink_particle_mass[-1][sink_i]*solar_mass), coord_system="data", text_args={"fontsize": 5, "color":'blue'})

                for p in [p1,p2]:
                    # Ensure the axes and colorbar limits match for all plots
                    plot = p.plots[field_name]
                    plot.figure = fig
                    if p==p1:
                        plot.axes = grid[i].axes
                    else:
                        plot.axes = grid[i+3].axes
                    plot.cax = grid.cbar_axes[i]

                    # Actually redraws the plot.
                    p.render()

                    # Modify the axes properties **after** p.render() so that they
                    # are not overwritten.
                i+=1
            
            plt.savefig("plots_to_sort/SinkVsNonSink_%s_%s/Core_%03d/DD%04d.png"%(field_name[1],which_cores,core_id,this_frame))
            plt.close()
    
    if set_clim==True:
        SinkVsNonSinkPlotter(trackname, sink_trackname, field_name = field_name, which_cores = which_cores, suppress_annotation = suppress_annotation, set_clim = False, clim_arr = [min_val,max_val])
        return True
    for core_id in main_core_list:
        animate_plots.animator("plots_to_sort/SinkVsNonSink_%s_%s/Core_%03d"%(field_name[1],which_cores,core_id),"Core_%03d"%(core_id))


SinkVsNonSinkPlotter('nb101', 'nb102', which_cores = 'densest', set_clim = False, clim_arr = [1e-3,1e4])
SinkVsNonSinkPlotter('nb101', 'nb102', field_name=("gas","magnetic_field_strength"), which_cores = 'densest', set_clim = False, clim_arr = [1e-2,1e5])
SinkVsNonSinkPlotter('nb101', 'nb102', field_name=("gas","velocity_magnitude"), which_cores = 'densest', set_clim = False, clim_arr = [1e-1,1e3])
SinkVsNonSinkPlotter('nb101', 'nb102', which_cores = 'alone', set_clim = False, clim_arr = [1e-3,1e4])
SinkVsNonSinkPlotter('nb101', 'nb102', field_name=("gas","magnetic_field_strength"), which_cores = 'alone', set_clim = False, clim_arr = [1e-2,1e5])
SinkVsNonSinkPlotter('nb101', 'nb102', field_name=("gas","velocity_magnitude"), which_cores = 'alone', set_clim = False, clim_arr = [1e-1,1e3])