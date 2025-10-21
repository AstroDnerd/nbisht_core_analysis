#import modules
import yt
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import warnings
import random
import sys
sys.path.append(os.path.abspath("/home/x-nbisht1/scripts"))
from animate_plots import animator
warnings.filterwarnings("ignore")
import json
from matplotlib.ticker import FuncFormatter, MultipleLocator
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid
#define constants

solar_mass = 5900*1.989e+33

def make_dir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def allsimplot(plot_dir, PlotType, path_to_data, dataset_to_plot, number_of_frames, start_frame = 0, cbar_range = None):
    #Define Arrays storing Sink Particle Data
    time = []
    #import data
    if not os.path.exists(path_to_output_plots+plot_dir+PlotType):
        os.makedirs(path_to_output_plots+plot_dir+PlotType)

    simnames = ['d03_Ms1.0', 'd07_Ms2.0', 'd11_Ms3.0', 'd15_Ms4.0', 'd19_Ms5.0', 'd23_Ms6.0', 'd27_Ms7.0', 'd31_Ms8.0']
    ddnames = ['DD0030', 'DD0050', 'DD0070', 'DD0090']
    for frame_number in range(start_frame,number_of_frames):
        frame_str = str(frame_number).zfill(4)
        fig = plt.figure(figsize=(5,8))
        # See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
        grid = AxesGrid(
            fig,
            (0.085, 0.085, 0.83, 0.83),
            nrows_ncols=(8, 4),
            axes_pad=0.0,
            label_mode="L",
            share_all=False,
            cbar_location="right",
            cbar_mode="single",
            cbar_size="3%",
            cbar_pad="0%",
            aspect=True,
        )
        for i,sim in enumerate(simnames):
            for j,dd in enumerate(ddnames):
                posn = i*4 + j
                simpath = path_to_data+sim+"_Ma0.0_512_"+dd+"/TT"+frame_str+"/time"+frame_str
                if os.path.exists(simpath):
                    print("Loading ", simpath)
                    data = yt.load(simpath)
                    plot = yt.ProjectionPlot(data, 'z', dataset_to_plot, fontsize = 10)
                    if i==0 and j==0:
                        plot.annotate_timestamp(corner="upper_left", draw_inset_box=True)
                    if i==0 and j==3:
                        plot.annotate_scale(corner="upper_right")
                    if cbar_range != None:
                        plot.set_zlim(dataset_to_plot[1],cbar_range[0],cbar_range[1])
                    ploty = plot.plots[dataset_to_plot]
                    ploty.figure = fig
                    ploty.axes = grid[posn].axes

                    ploty.cax = grid.cbar_axes[posn]

                    # Actually redraws the plot.
                    plot.render()

                    # Modify the axes properties **after** p.render() so that they
                    # are not overwritten.
                # Only show tick labels for bottom-left panel (i=7, j=0)
                if not (i == 7 and j == 0):
                    grid[posn].axes.set_xticklabels([])
                    grid[posn].axes.set_yticklabels([])
                if i==0:
                    grid[posn].axes.set_title(dd, fontsize=15)
                if j==0:
                    grid[posn].axes.set_ylabel(sim[4:], fontsize=15)
        
        
        plt.savefig(path_to_output_plots+plot_dir+PlotType+"/DD"+frame_str+".png", dpi=300)
        plt.close()

    animator(path_to_output_plots+plot_dir+PlotType,PlotType, duration = 0.5 )


path_to_data_512_nonsim = "/home/x-nbisht1/scratch/projects/512/NonsinkSimSuite/"
path_to_output_plots = "/home/x-nbisht1/projects/results/512/NonsinkSimSuite/"

reg_plot = 1
if reg_plot ==1:
    num_of_frames = 100
    allsimplot('Plots_allsim/','DensityPlot_all', path_to_data_512_nonsim, ("enzo", "Density"), num_of_frames, start_frame=0, cbar_range=None)


