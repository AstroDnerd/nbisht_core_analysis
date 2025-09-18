
import matplotlib
matplotlib.use('Agg')
import yt
import math
import matplotlib.pyplot as plt
import numpy as np
import pdb
nar = np.array

from importlib import reload

import looper
reload(looper)

bad_id = 1499257
if 'this_looper' not in dir():
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
    this_looper = looper.core_looper(directory= directory,
                                     sim_name = 'u05',
                                     out_prefix = 'test',
                                     target_frame = 125,
                                     frame_list = [50,55,60],#[0,1,2]+list(range(10,130,5))+[125],
                                     core_list = [31],
                                     fields_from_grid=['x','y','z']
                                  )
    this_looper.get_target_indices(h5_name='u05_0125_peaklist.h5',
                                     bad_particle_list='bad_particles.h5')
    this_looper.target_indices[31]=this_looper.target_indices[31][::100]
    if bad_id not in this_looper.target_indices[31]:
        this_looper.target_indices[31][0]=bad_id
    #this_looper.target_indices={31:np.array([1499257])}
    #this_looper.target_indices={31:np.array([1499257,  145997,  207822, 1368564])}#, 1685267, 1741447, 2064187])}
    this_looper.get_tracks()

import pyximport; pyximport.install()
import particle_ops
import particle_grid_mask
#examine get_current_mask

if 'data_region' not in dir():
    these_pids=np.array([bad_id]) #self.target_indices.astype('int64')
    #data_region = self.get_region(self.frame)
    data_region=this_looper.snaps[60][31].ds.all_data()
    mask_to_get=np.zeros(these_pids.shape,dtype='int32')
    my_indices = data_region['particle_index'].astype('int64')
bad_index = np.where(my_indices==bad_id)
if 'gtots' not in dir():
    gtots = this_looper.snaps[60][31].ds.index.grids[66]
    snap = this_looper.snaps[60][31]


if 'getmask_pos' not in dir():
    getmask_pos = data_region['particle_position'][bad_index]
    getmask_vel = data_region['particle_velocity'][bad_index]
    grid_index=75
    print( 'got the right one: ',gtots['particle_index'][grid_index]==bad_id)
    grid_pos = gtots['particle_position'][grid_index]
    print( getmask_pos-grid_pos)

if 'found_any' not in dir():
    found_any, mask = particle_ops.mask_particles(these_pids,my_indices,mask_to_get)
    data_region['particle_index'][mask.astype('bool')] == bad_id
if 1:
    index_3= np.where( snap.ind == bad_id)
    print("pos: mask",getmask_pos)
    print("pos: snap",snap.pos[index_3])
    print("den: fv", snap.field_values['density'][index_3])

if 1:
    thtr = this_looper.tr
    track_id = np.where(thtr.particle_ids==1499257)
    tr_density = thtr.c(31,'density')
    this_density = tr_density[track_id,-1]
    #print("The density for that particle",this_density)
    print("The density for that particle",thtr.p([bad_id],'density'))

if 1:
    #so far everyone agrees.  What's the problem again?
    if 0:
        print(this_looper.snaps[60][31].pos)
        print(getmask_pos)
        print(this_looper.snaps[60][31].vel)
        print(getmask_vel)
    fv = this_looper.snaps[60][31].field_values
    ind = np.where( snap.ind == bad_id)
    if 0:
        print("%0.16f %0.16f %0.16f"%(fv['x'][ind],fv['y'][ind],fv['z'][ind]))
        print("%0.16f %0.16f %0.16f"%tuple(list(getmask_pos.v[0])) )
    zones = np.array( [ (fv['x'][ind]-gtots.LeftEdge[0].v)/gtots.dds[0],
                        (fv['y'][ind]-gtots.LeftEdge[1].v)/gtots.dds[1],
                        (fv['z'][ind]-gtots.LeftEdge[2].v)/gtots.dds[2]])
    bzones = (getmask_pos.v-gtots.LeftEdge.v)/gtots.dds.v
    tr_x = thtr.p([bad_id],'x')
    tr_y = thtr.p([bad_id],'y')
    tr_z = thtr.p([bad_id],'z')
    print(tr_x-fv['x'][ind])
    print(tr_y-fv['y'][ind])
    print(tr_z-fv['z'][ind])

    
#found_any, mask = particle_ops.mask_particles(these_pids,my_indices,mask_to_get)
#self.mask = mask
