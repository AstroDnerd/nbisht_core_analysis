from starter1 import *
import matplotlib
matplotlib.use('Agg')
import yt
import math
import matplotlib.pyplot as plt
import numpy as np
import pdb
nar = np.array
fptr = open('n_particles.txt','r')
lines=fptr.readlines()
fptr.close()
parts = np.zeros([len(lines),2])
for n,line in enumerate(lines):
    parts[n] = np.array(line.split(),dtype='int')
all_nonzero = parts[:,0][ parts[:,1] >0]
from importlib import reload

import looper
reload(looper)
import trackage
reload(trackage)
import tracks_read_write as trw
reload(trw)
#many1: frame_list = list(range(10,130,10)) + [125]
#many1: core_list =  [ 120, 10, 308, 41, 44, 110, 1, 65], 
#frame 120 core 96
core_list = all_nonzero.astype('int')
if 'this_looper' not in dir():
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
    this_looper = looper.core_looper(directory= directory,
                                     sim_name = 'u05',
                                     out_prefix = 'test',
                                     target_frame = 125,
                                     frame_list = [1,100,120],#[0,1,2]+list(range(10,130,10))+[125],
                                     core_list =  [],# core_list,
                                     fields_from_grid=['x','y','z','velocity_magnitude','magnetic_field_strength',
                                                      'velocity_divergence']
                                  )
for nfile, fname in enumerate(glob.glob('track_three_to_test*.h5')):
    trw.load_loop(this_looper,fname)
