nonsink_trackname = 'nb101'
sink_trackname = 'nb102'
#import modules
import yt
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import warnings
import random
import sys
from animate_plots import animator
warnings.filterwarnings("ignore")
import json


#define paths sinkoff
path_to_data_128 = '/data/cb1/nbisht/anvil_scratch/projects/128/B2/'


path_to_output_plots = "plots_to_sort/"


#define constants


solar_mass = 1.989e+33

def make_dir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def regular_plot(plot_dir, PlotType, path_to_data, coordinate, dataset_to_plot, number_of_frames,cbar_range = None):
    #import data
    if not os.path.exists(path_to_output_plots+plot_dir+PlotType):
        os.makedirs(path_to_output_plots+plot_dir+PlotType)

    for frame_number in range(0,number_of_frames):
        frame_str = str(frame_number).zfill(4)
        data = yt.load(path_to_data+"DD"+frame_str+"/data"+frame_str)
        all_data = data.all_data()
        print(data.field_list)

        plot = yt.ProjectionPlot(data, coordinate, dataset_to_plot, fontsize = 10)
        plot.annotate_timestamp(corner="upper_left", redshift=True, draw_inset_box=True)
        plot.annotate_scale(corner="upper_right")
        plot.annotate_particles((1,'cm'), p_size=0.1)
        if cbar_range != None:
            plot.set_zlim(dataset_to_plot[1],cbar_range[0],cbar_range[1])
        plot.save(path_to_output_plots+plot_dir+PlotType+"/DD"+frame_str+".png")

    animator(path_to_output_plots+plot_dir+PlotType,PlotType)


reg_plot = 1
if reg_plot ==1:
    num_of_frames = 136
    regular_plot('NonSink_Density/','Particle_xy', path_to_data_128, 'z', ("enzo", "Density"), num_of_frames,cbar_range=None)
    regular_plot('NonSink_Density/','Particle_yz', path_to_data_128, 'x', ("enzo", "Density"), num_of_frames,cbar_range=None)
