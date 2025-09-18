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

def SinkVsNonSink_BRhoPlotter(trackname, sink_trackname, framelist):
    make_dir("plots_to_sort/B_Rho_phasediag")
    make_dir("plots_to_sort/B_Rho_phasediag_diff")
    for this_frame in framelist:
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
            

        plt.savefig("plots_to_sort/B_Rho_phasediag/DD%04d.png"%(this_frame))
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
        plt.savefig("plots_to_sort/B_Rho_phasediag_diff/DD%04d.png"%(this_frame))
        plt.close()
    
    animate_plots.animator("plots_to_sort/B_Rho_phasediag","B_Rho_phasediag")
    animate_plots.animator("plots_to_sort/B_Rho_phasediag_diff","B_Rho_phasediag_diff")



#SinkVsNonSink_BRhoPlotter('nb101', 'nb102', np.arange(0,126))