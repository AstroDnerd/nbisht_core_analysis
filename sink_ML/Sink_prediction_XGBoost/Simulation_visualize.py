#import modules
import yt
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import warnings
import random
import sys

warnings.filterwarnings("ignore")
import json


#define paths
SIMDIR_nb101  = '/data/cb1/nbisht/anvil_scratch/projects/128/B2/'
SIMDIR_u502  = '/data/cb1/Projects/P19_CoreSimulations/u202-Beta2/GravPotential/'

path_to_output_plots = './plots_to_sort/SimViz/'

#define constants


solar_mass = 1.989e+33

def make_dir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def regular_plot(plot_dir, PlotType, path_to_data, coordinate, dataset_to_plot, number_of_frames,cbar_range = None, zoom_level = 1):
    #import data
    if not os.path.exists(path_to_output_plots+plot_dir+PlotType):
        os.makedirs(path_to_output_plots+plot_dir+PlotType)

    for frame_number in range(0,number_of_frames):
        frame_str = str(frame_number).zfill(4)
        data = yt.load(path_to_data+"DD"+frame_str+"/data"+frame_str)
        all_data = data.all_data()
        print(data.field_list)

        plot = yt.ProjectionPlot(data, coordinate, dataset_to_plot, fontsize = 10)
        plot.set_cmap(dataset_to_plot, "dusk")
        plot.set_background_color(dataset_to_plot, color="black")
        plot.zoom(zoom_level)
        plot.annotate_timestamp(corner="upper_left", redshift=True, draw_inset_box=True)
        plot.annotate_scale(corner="upper_right")
        if cbar_range != None:
            plot.set_zlim(dataset_to_plot[1],cbar_range[0],cbar_range[1])
        plot.save(path_to_output_plots+plot_dir+PlotType+"/DD"+frame_str+".png")


def volume_renderer(plot_dir, PlotType, path_to_data, coordinate, dataset_to_plot, number_of_frames,cbar_range = None, zoom_level = 1):
    #import data
    if not os.path.exists(path_to_output_plots+plot_dir+PlotType):
        os.makedirs(path_to_output_plots+plot_dir+PlotType+"/transfer_function")
        os.makedirs(path_to_output_plots+plot_dir+PlotType+"/volume_rendering")

    for frame_number in range(0,number_of_frames):
        frame_str = str(frame_number).zfill(4)
        data = yt.load(path_to_data+"DD"+frame_str+"/data"+frame_str)
        all_data = data.all_data()
        print(data.field_list)

        sc = yt.create_scene(data, lens_type="perspective")

        source = sc[0]

        source.set_field(dataset_to_plot)
        source.set_log(True)

        bounds = (1e-4, 1e4)

        # Since this rendering is done in log space, the transfer function needs
        # to be specified in log space.
        tf = yt.ColorTransferFunction(np.log10(bounds))

        tf.add_layers(5, colormap="dusk")

        source.tfh.tf = tf
        source.tfh.bounds = bounds

        source.tfh.plot(path_to_output_plots+plot_dir+PlotType+"/transfer_function/DD"+frame_str+".png", profile_field=("gas", "density"))

        sc.save(path_to_output_plots+plot_dir+PlotType+"/volume_rendering/DD"+frame_str+".png")#, sigma_clip=6)



reg_plot = 0
if reg_plot ==1:
    num_of_frames = 115
    regular_plot('DensityPlots/','DensityPlot_xy_nb101', SIMDIR_nb101, 'z', ("enzo", "Density"), num_of_frames,cbar_range=None,zoom_level = 1)
    regular_plot('DensityPlots/','DensityPlot_xy_u502', SIMDIR_u502, 'z', ("enzo", "Density"), num_of_frames,cbar_range=None,zoom_level = 1)

vol_rend = 1
if vol_rend ==1:
    num_of_frames = 115
    volume_renderer('DensityPlots/','DensityPlot_3D_render_nb101', SIMDIR_nb101, 'z', ("enzo", "Density"), num_of_frames,cbar_range=None,zoom_level = 1)