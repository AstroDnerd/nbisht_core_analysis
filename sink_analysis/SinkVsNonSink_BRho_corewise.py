import warnings
warnings.filterwarnings("ignore")
import h5py
import numpy as np
import random
import pandas
from starter2 import *
import json
from mpl_toolkits.axes_grid1 import AxesGrid
import animate_plots
import track_loader as TL

def SinkVsNonSink_BRhoPlotter_corewise(trackname, sink_trackname, which_cores = 'densest'):
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
        this_track = track_info.tracks[trackname]
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
    for this_frame_index in range(len(tracklist)):
        this_frame = tracklist[this_frame_index]
        this_track = track_info.tracks[trackname]
        this_simname = this_track.sim_directory

        ds = yt.load("%s/DD%04d/data%04d"%(this_simname,this_frame,this_frame))
        all_data = ds.all_data()

        #load sink simulation
        sink_track = track_info.tracks[sink_trackname]
        sink_simname = sink_track.sim_directory
        sink_frame = this_frame
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

        weight_arr = ["cell_volume", "cell_mass", "velocity_magnitude"]
        weight_lim = [[1e-10,1e-6],[1e-11,1e-1],[1e-1,1e+2]]
        for core_id in main_core_list:
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
            i = 0
            j=0
            z_data_arr = []
            make_dir("plots_to_sort/B_Rho_phasediag_corewise_%s/Core_%03d"%(which_cores,core_id))
            make_dir("plots_to_sort/B_Rho_phasediag_diff_corewise_%s/Core_%03d"%(which_cores,core_id))
            ms = mini_scrubbers[core_id]
            print('Frame %04d, core %03d'%(this_frame,core_id))
            ms.particle_pos(core_id)
            core_center = ms.mean_center[:,this_frame_index]
            core_rms = ms.rms[this_frame_index]
            core_radius = core_rms+0.01
            ds_sphere = ds.sphere(core_center, core_radius)
            ds_sink_sphere = ds_sink.sphere(core_center, core_radius)
            for Brho_weight in weight_arr:   
                p1 = yt.PhasePlot(ds_sphere,("gas","density"),("gas","magnetic_field_strength"), ("gas",Brho_weight), fontsize = 10)
                p2 = yt.PhasePlot(ds_sink_sphere,("gas","density"),("gas","magnetic_field_strength"), ("gas",Brho_weight), fontsize = 10)
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
                    z_data_arr.append(z_data.T)
                    plot.axes.contour(x_cen,y_cen,z_data)
                    if p==p2:
                        norm = mpl.colors.Normalize(vmin=weight_lim[j][0], vmax=weight_lim[j][1])
                        for sink_i in range(len(sink_particle_mass[-1])):
                            sink_pos = [float(sink_particle_xpos[-1][sink_i]), float(sink_particle_ypos[-1][sink_i]),float(sink_particle_zpos[-1][sink_i])]
                            new_x = float(sink_particle_xpos[-1][sink_i]) - core_center[0]
                            new_y = float(sink_particle_ypos[-1][sink_i]) - core_center[1]
                            new_z = float(sink_particle_zpos[-1][sink_i]) - core_center[2]
                            if (new_x<core_radius and new_x>-core_radius) and (new_y<core_radius and new_y>-core_radius) and (new_z<core_radius and new_z>-core_radius):
                                point_obj =  ds_sink.point(sink_pos)
                                sink_cell_volume = point_obj['gas', 'cell_volume']
                                sink_cell_density = sink_particle_mass[-1][sink_i]/float(sink_cell_volume)
                                sink_cell_magfield = float(point_obj['gas', 'magnetic_field_strength'])
                                color_variable_val = float(point_obj['gas', Brho_weight])
                                c_sink = plot.cb.cmap(norm(color_variable_val))
                                plot.axes.scatter(sink_cell_density,sink_cell_magfield, s = 5, c = c_sink)

                    i+=1
                j+=1
                

            plt.savefig("plots_to_sort/B_Rho_phasediag_corewise_%s/Core_%03d/DD%04d.png"%(which_cores,core_id,this_frame))
            plt.close()

            fig = plt.figure()
            grid = AxesGrid(
                fig,
                (0.085, 0.085, 0.9, 0.9),
                nrows_ncols=(1, 3),
                axes_pad=0.5,
                label_mode="L",
                share_all=True,
                aspect=True,
            )
            im = grid[0].axes.imshow(z_data_arr[0]-z_data_arr[1], cmap = 'PRGn')
            fig.colorbar(im, ax=grid[0].axes, location='top', anchor=(0, -1.2), shrink=0.3)
            im = grid[1].axes.imshow(z_data_arr[2]-z_data_arr[3], cmap = 'PRGn')
            fig.colorbar(im, ax=grid[1].axes, location='bottom', anchor=(0.5, 1.7), shrink=0.3)
            im = grid[2].axes.imshow(z_data_arr[4]-z_data_arr[5], cmap = 'PRGn')
            fig.colorbar(im, ax=grid[2].axes, location='top', anchor=(1, -1.2), shrink=0.3)
            plt.savefig("plots_to_sort/B_Rho_phasediag_diff_corewise_%s/Core_%03d/DD%04d.png"%(which_cores,core_id,this_frame))
            plt.close()
    
    for core_id in main_core_list:
        animate_plots.animator("plots_to_sort/B_Rho_phasediag_corewise_%s/Core_%03d"%(which_cores,core_id),"Core_%03d"%(core_id))
        animate_plots.animator("plots_to_sort/B_Rho_phasediag_diff_corewise_%s/Core_%03d"%(which_cores,core_id),"Core_%03d"%(core_id))



SinkVsNonSink_BRhoPlotter_corewise('nb101', 'nb102', which_cores = 'alone')