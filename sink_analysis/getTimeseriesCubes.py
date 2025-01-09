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
import scipy.fft as fft
from animate_plots import animator
import pdb
from starter2 import *
import track_loader as TL


nonsink_trackname = 'nb101'
sink_trackname = 'nb102'

#define paths sinkoff
path_to_output_plots = "./plots_to_sort/TimeCorrelation/"

def make_dir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


def getTimeseriesCubes(trackname, fieldname = 'density', target_frames = None):
    this_track = track_info.tracks[trackname]
    import time
    start_time = time.time()
    if target_frames ==None:
        target_frames = this_track.frame_list
    
    df_name = this_track.sim_directory+'/datasets/TimeseriesCubes_Density.json'
    if os.path.isfile(df_name):
        print("File exists!")
        return 1

    TSCube_array = []
    time_arr = []
    for framenum in target_frames:
        #Open these frames and get particle data
        ds_begin = yt.load("%s/DD%04d/data%04d"%(this_track.sim_directory,framenum,framenum))
        tcur = float(ds_begin.current_time)
        if time_arr!=[]:
            if tcur==time_arr[-1]:
                continue
        all_data_begin = ds_begin.all_data()
        level = 0
        cg = ds_begin.covering_grid(level, left_edge=ds_begin.domain_left_edge, dims=ds_begin.domain_dimensions * 2**level, fields=[fieldname])
        TSCube_array.append(cg[fieldname].tolist())
        time_arr.append(float(ds_begin.current_time))
    
    TSCube_dic = {"Cube":TSCube_array, "Time": time_arr}
    print(time_arr)
    infile = open(df_name, 'w')
    json.dump(TSCube_dic, infile)
    infile.close()
    diff = time.time() - start_time
    print('Time taken by looper:',diff)

getTimeseriesCubes(nonsink_trackname, fieldname = 'density')


