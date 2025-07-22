#import modules
import yt
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import warnings
import random
import sys
sys.path.insert(0, '/home/nbisht/nikhilb_home/scripts/nbisht_core_analysis/')
warnings.filterwarnings("ignore")
import json
import scipy.fft as fft
import pdb
from starter2 import *
import track_loader as TL

'''
ENZO FIELD LIST
[('all', 'particle_index'),
 ('all', 'particle_mass'),
 ('all', 'particle_position_x'),
 ('all', 'particle_position_y'),
 ('all', 'particle_position_z'),
 ('all', 'particle_type'),
 ('all', 'particle_velocity_x'),
 ('all', 'particle_velocity_y'),
 ('all', 'particle_velocity_z'),
 ('enzo', 'Bx'),
 ('enzo', 'BxF'),
 ('enzo', 'By'),
 ('enzo', 'ByF'),
 ('enzo', 'Bz'),
 ('enzo', 'BzF'),
 ('enzo', 'Dark_Matter_Density'),
 ('enzo', 'Density'),
 ('enzo', 'DivB'),
 ('enzo', 'Ex'),
 ('enzo', 'Ey'),
 ('enzo', 'Ez'),
 ('enzo', 'PotentialField'),
 ('enzo', 'Temperature'),
 ('enzo', 'x-velocity'),
 ('enzo', 'y-velocity'),
 ('enzo', 'z-velocity'),
 ('io', 'particle_index'),
 ('io', 'particle_mass'),
 ('io', 'particle_position_x'),
 ('io', 'particle_position_y'),
 ('io', 'particle_position_z'),
 ('io', 'particle_type'),
 ('io', 'particle_velocity_x'),
 ('io', 'particle_velocity_y'),
 ('io', 'particle_velocity_z'),
 ('nbody', 'particle_index'),
 ('nbody', 'particle_mass'),
 ('nbody', 'particle_position_x'),
 ('nbody', 'particle_position_y'),
 ('nbody', 'particle_position_z'),
 ('nbody', 'particle_type'),
 ('nbody', 'particle_velocity_x'),
 ('nbody', 'particle_velocity_y'),
 ('nbody', 'particle_velocity_z')]
'''
nonsink_trackname = 'nb101'
sink_trackname = 'nb102'

#define paths sinkoff
path_to_output_plots = "./plots_to_sort/TimeCorrelation/"

def make_dir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


def getTimeseriesCubes(trackname, fieldname = 'density', target_frames = None, savename = None):
    this_track = track_info.tracks[trackname]
    import time
    start_time = time.time()
    if target_frames ==None:
        target_frames = this_track.frame_list
    
    df_name = track_info.tracks['nb101'].sim_directory+'/datasets/'+trackname+'_TimeseriesCubes_'+savename+'.npy'
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
    infile = open(df_name, 'wb')
    np.save(infile, TSCube_array)
    np.save(infile, time_arr)
    infile.close()
    diff = time.time() - start_time
    print('Time taken by looper:',diff)


getTimeseriesCubes('u503', fieldname = 'Density', savename = 'Density')
getTimeseriesCubes('u503', fieldname = 'x-velocity', savename = 'velocity_x')
getTimeseriesCubes('u503', fieldname = 'y-velocity', savename = 'velocity_y')
getTimeseriesCubes('u503', fieldname = 'z-velocity', savename = 'velocity_z')


