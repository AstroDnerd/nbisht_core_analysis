
from starter2 import *
import xtra_energy
import data_locations as dl
reload(dl)
if 'this_simname' not in dir():
    this_simname = 'u202'

frame_list=[0]
fields=['density']

# to get this range, do h5ls on the .h5
core_list=range(408)
if 1:
    this_looper = looper.core_looper(directory= dl.sims[this_simname],
                                     sim_name = this_simname,
                                     out_prefix = this_simname,
                                     target_frame = dl.target_frames[this_simname],
                                     frame_list = frame_list,
                                     core_list =  core_list,
                                     fields_from_grid=fields
                                  )
    this_looper.get_target_indices(h5_name=dl.peak_list[this_simname])
                                     #bad_particle_list=dl.bad_particles[this_simname])

def count_particles(this_looper,fname='n_particles.txt'):
    fptr = open(fname,'w')
    for core in this_looper.target_indices:
        fptr.write("%d %d\n"%(core, len(this_looper.target_indices[core])))
    fptr.close()
count_particles(this_looper,'%s_n_particles.txt'%this_simname)
