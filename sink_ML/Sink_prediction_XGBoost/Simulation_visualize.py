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

from mpl_toolkits.axes_grid1 import AxesGrid
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


def paper_regular_plot(plot_dir, PlotName, path_to_data, coordinate, dataset_to_plot, number_of_frames):
    fig = plt.figure(figsize=(15, 4))

    grid = AxesGrid(
        fig,
        (0.075, 0.075, 0.85, 0.85),
        nrows_ncols=(1, 4),
        axes_pad=0.05,
        label_mode="L",
        share_all=True,
        cbar_location="right",
        cbar_mode="single",
        cbar_size="3%",
        cbar_pad="0%",
    )
    i=0
    for frame_number in [0,25,75,num_of_frames-1]:
        frame_str = str(frame_number).zfill(4)
        data = yt.load(path_to_data+"DD"+frame_str+"/data"+frame_str)
        all_data = data.all_data()
        print(data.field_list)

        p = yt.ProjectionPlot(data, coordinate, dataset_to_plot, fontsize = 10)
        p.set_cmap(dataset_to_plot, "CMRmap")
        p.set_background_color(dataset_to_plot, color="black")
        p.annotate_timestamp(corner="upper_left")
        p.annotate_scale(corner="upper_right")
        p.set_zlim(dataset_to_plot[1],1e-2,1e4)
        plot = p.plots[dataset_to_plot]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
        p.render()
        i+=1
    plt.savefig(path_to_output_plots+plot_dir+PlotName+".png", dpi=300, bbox_inches="tight")



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

vol_rend = 0
if vol_rend ==1:
    num_of_frames = 115
    volume_renderer('DensityPlots/','DensityPlot_3D_render_nb101', SIMDIR_nb101, 'z', ("enzo", "Density"), num_of_frames,cbar_range=None,zoom_level = 1)

paper_reg_plot = 1
if paper_reg_plot ==1:
    num_of_frames = 115
    #paper_regular_plot('DensityPlots/','DensityPlot_xy_nb101', SIMDIR_nb101, 'z', ("enzo", "Density"), num_of_frames)
    paper_regular_plot('DensityPlots/','DensityPlot_xy_u502', SIMDIR_u502, 'z', ("enzo", "Density"), num_of_frames)