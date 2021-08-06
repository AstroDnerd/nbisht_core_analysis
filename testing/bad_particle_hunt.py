
import matplotlib
matplotlib.use('Agg')
import yt
import math
import matplotlib.pyplot as plt
import numpy as np
import pdb
nar = np.array
#fptr = open('n_particles.txt','r')
#lines=fptr.readlines()
#fptr.close()
#parts = np.zeros([len(lines),2])
#for n,line in enumerate(lines):
#    parts[n] = np.array(line.split(),dtype='int')
#all_nonzero = parts[:,0][ parts[:,1] >0]
#from importlib import reload


def check_particles(ds):
    bad_index=[]
    for grid in ds.index.grids:
        pos = grid['particle_position']
        for i in [0,1,2]:
            check=pos[:,i] > grid.RightEdge[i]
            check=np.logical_or(check,pos[:,i] < grid.LeftEdge[i])
            if check.any():
                locations=np.where(check)
                bad_index+= list( grid['particle_index'][check].v)
    return bad_index
#g=check_particles(this_looper.ds_list[120])
bad_list = []
#frame_list = list(range(122))+[125]
this_simname = 'u201'
frame_list = range( dl.target_frames[this_simname]+1)
directory= dl.sims[this_simname]

for frame in frame_list:
    ds = yt.load("%s/DD%04d/data%04d"%(directory,frame,frame))
    this_bad_list = check_particles(ds)
    print("Frame %d len bad %d"%(frame,len(this_bad_list)))
    bad_list += this_bad_list
    bad_list = list(np.unique(bad_list))

