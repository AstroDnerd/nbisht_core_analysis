import matplotlib
matplotlib.use('Agg')
import math
import matplotlib.pyplot as plt
import numpy as np
nar = np.array
import time
from importlib import reload
import h5py
import copy 
import pdb
from importlib import reload
from collections import defaultdict
import weakref
if 0:
    import loop_tools
    reload(loop_tools)
"""
trackage.
the tool package to get history for individual particles
"""


class track_manager():
    """Manages tracks.
    track_manager[field] returns a list of [particle, time]
    so 
    for i in track_manager[field]:
        plt.plot(track_manager.time,i)
    plots the field vs. time.
    self.ingest(snapshot) populates the master list
    """
    def __init__(self,my_loop=None):
        self.my_loop = None
        if my_loop is not None:
            self.my_loop = weakref.proxy(my_loop)
        self.particle_ids=nar([],dtype='int64')
        self.core_ids = nar([],dtype='int64')
        self.frames = nar([],dtype='int64')
        self.times = nar([],dtype='float64')
        self.V_rad = nar([],dtype='float64')
        self.V_rel = nar([],dtype='float64')
        self.track_dict={}
        self.shape=(0,0) #particles, frames
    def sort_time(self):
        """fix up the time ordering of the members"""
        asort =  np.argsort(self.times)
        if (asort != sorted(asort)).any():
            for k in self.track_dict:
                self.track_dict[k]=self.track_dict[k][:,asort]
            self.times = self.times[asort]
            self.frames= self.frames[asort]

    def write(self,fname):
        fptr = h5py.File(fname,'w')
        try:
            fptr['particle_ids'] = self.particle_ids
            fptr['core_ids'] = self.core_ids
            fptr['frames'] = self.frames
            fptr['times'] = self.times
            fptr['r_velocity'] = self.V_rad
            fptr['full_velocity']=self.V_rel
            for k in self.track_dict:
                fptr[k] = self.track_dict[k]
        except:
            raise
        finally:
            fptr.close()
    def read(self,fname):
        fptr = h5py.File(fname,'r')
        try:
            for key in fptr:
                if str(key) in ['particle_ids', 'core_ids', 'frames', 'times','r_velocity','full_velocity']:
                    self.__dict__[key]=fptr[key][:]
                else:
                    self.track_dict[key] = fptr[key][:] 
        except:
            raise
        finally:
            fptr.close()
    def merge(self,fname):
        fptr = h5py.File(fname,'r')
        temp_dict = {}
        temp_track_dict = {}
        try:
            for key in fptr:
                if str(key) in ['particle_ids', 'core_ids', 'frames', 'times']:
                    temp_dict[key]=fptr[key][:]
                else:
                    temp_track_dict[key] = fptr[key][:] 
        except:
            raise
        finally:
            fptr.close()
        if len(temp_dict['particle_ids']) != len(self.particle_ids):
            print("Error with merging: wrong number of particles")
        elif np.abs(temp_dict['particle_ids'] - self.particle_ids).sum() > 0:
            print("error with merging: wrong particles")
        elif np.abs( temp_dict['core_ids']-self.core_ids).sum() > 0:
            print("error with merging: core ids")
        else:
            frames = np.concatenate([self.frames,temp_dict['frames']])
            args= np.argsort(frames)
            self.frames=frames[args]
            self.times=np.concatenate([self.times,temp_dict['times']])[args]
            for key in self.track_dict:
                oot = np.concatenate( [self.track_dict[key], temp_track_dict[key]],axis=1)[:,args]
                self.track_dict[key]=oot

        return {'td':temp_dict,'ttd':temp_track_dict}

    def ingest(self,snapshot):
        #pdb.set_trace()
        particle_ids = copy.copy(self.my_loop.target_indices[snapshot.core_id])
        if snapshot.core_id not in self.core_ids:
            #this might not be the best place for the parent step.
            core_ids = np.ones_like(particle_ids) * snapshot.core_id
            self.core_ids = np.append(self.core_ids, core_ids)
            self.particle_ids = np.append(self.particle_ids, particle_ids)
        particle_start = np.where(self.particle_ids==particle_ids[0])[0][0]
        particle_end=particle_start+particle_ids.size

        if snapshot.frame not in self.frames:
            self.frames=np.append(self.frames,snapshot.frame)
            self.times=np.append(self.times,snapshot.time) #ds['InitialTime'])
            self.V_rad = np.append(self.V_rad,snapshot.V_radial)

        frame_id = np.where(self.frames == snapshot.frame)[0][0]

        for field in snapshot.field_values:
            current_shape = self[field].shape
            new_shape = [self.particle_ids.size,
                         self.frames.size]
            temp_frame = np.zeros(new_shape,dtype='float64')
            old_slice = (slice(None,current_shape[0]),
                         slice(None,current_shape[1]))
            temp_frame[old_slice]= self[field]
            new_slice = (slice(particle_start,particle_end),
                         slice(frame_id,frame_id+1))
            nuggle=np.array(snapshot.field_values[field])
            nuggle.shape=(particle_ids.size,1)
            temp_frame[new_slice]=nuggle
            self[field]=temp_frame


    def p(self,particle_list,field):
        output = None
        for particle in particle_list:
            loc = np.where(self.particle_ids == particle)
            parts = self[field][loc,:]
            if output is None:
                output = parts
            else:
                output= np.append(output, parts,axis=0)
        return output
    def c(self,core_list,field):
        if type(core_list) is int:
            core_list = [core_list]
        nf = self.frames.size
        output = None
        for core in core_list:
            loc = self.core_ids == core
            core_values = self[field][loc,:]
            if output is None:
                output = core_values
            else:
                output = np.append(output,core_values,axis=0)
        return output
    def __getitem__(self,item):
        output = None
        if item in ['particles', 'particle_ids']:
            return self.particle_ids
        if item in ['core_ids','cores']:
            return self.core_ids
        if item not in self.track_dict:
            self.track_dict[item]=np.zeros(self.shape,dtype='float64')
        return self.track_dict[item]
    def __setitem__(self,item,value):
        self.track_dict[item]=value
    def setup_fields(self,field_list=None):
        for frame in self.my_loop.field_list:
            for core_id in self.my_loop.core_list:
                this_snapshot = looper.make_snapshot(frame,core_id)
                if field_list is None:
                    local_field_list = this_snapshot.field_values.keys()
                else:
                    local_field_list = field_list
                for field in local_field_list:
                    self[field].ingest(this_snapshot)


def shift_down(pos):
    #shift based on the time history: if a particle jumps more than half the box,
    #move it.
    sign = -1
    out = copy.copy(pos)

    mean_pos = np.median(out[:,-1])
    #print("mean_pos",mean_pos)
    distance_from_final =        out- mean_pos
    #ft = np.abs(distance_from_final[:,:-1]) > 0.5
    ft = np.abs(distance_from_final) > 0.5
    out[ft] -=  1*np.sign(distance_from_final[ft])



    #out[ out > point ] = out[ out > point]-1
    return out#,delta
        
class mini_scrubber():
    def __init__(self,trk,core_id,do_velocity=False):
        self.trk=trk
        self.scrub(core_id,do_velocity=do_velocity)
        self.axis=0
                
    def scrub(self,core_id, axis=0, do_velocity=False):
        self.raw_x = self.trk.c([core_id],'x')
        self.raw_y = self.trk.c([core_id],'y')
        self.raw_z = self.trk.c([core_id],'z')

        self.density = self.trk.c([core_id],'density')
        self.cell_volume = self.trk.c([core_id],'cell_volume')
        self.mass = self.density*self.cell_volume
        self.mass_total=self.mass.sum(axis=0)
        #this_x=raw_x
        #this_y=raw_y
        if 1:
            #do the shift
            self.this_x = shift_down(self.raw_x)
            self.this_y = shift_down(self.raw_y)
            self.this_z = shift_down(self.raw_z)
        else:
            #don't actuall shift
            self.this_x = self.raw_x+0
            self.this_y = self.raw_y+0
            self.this_z = self.raw_z+0
        self.mean_x = np.sum(self.this_x*self.mass,axis=0)/self.mass_total
        self.mean_y = np.sum(self.this_y*self.mass,axis=0)/self.mass_total
        self.mean_z = np.sum(self.this_z*self.mass,axis=0)/self.mass_total
        self.nparticles,self.ntimes=self.this_x.shape
        self.meanx2 = np.tile(self.mean_x,(self.raw_x.shape[0],1))
        self.meany2 = np.tile(self.mean_y,(self.raw_x.shape[0],1))
        self.meanz2 = np.tile(self.mean_z,(self.raw_z.shape[0],1))
        self.rx_rel=self.this_x-self.meanx2
        self.ry_rel=self.this_y-self.meany2
        self.rz_rel=self.this_z-self.meanz2
        self.r2 = self.rx_rel**2+self.ry_rel**2+self.rz_rel**2
        self.r=np.sqrt(self.r2)
        self.rmax = np.max(self.r,axis=0)
        self.max_track = np.where( self.r[:,0] == self.rmax[0])
        self.rmax_fat=np.tile(self.rmax,(self.raw_x.shape[0],1))
        self.rms = np.sqrt( np.mean(self.r2,axis=0))

        if do_velocity:
            self.raw_vx = self.trk.c([core_id],'velocity_x')
            self.raw_vy = self.trk.c([core_id],'velocity_y')
            self.raw_vz = self.trk.c([core_id],'velocity_z')
            self.mean_vx = np.sum(self.raw_vx*self.mass,axis=0)/self.mass_total
            self.mean_vy = np.sum(self.raw_vy*self.mass,axis=0)/self.mass_total
            self.mean_vz = np.sum(self.raw_vz*self.mass,axis=0)/self.mass_total
            self.sqr_vx = np.sum(self.raw_vx**2,axis=0)
            self.sqr_vy = np.sum(self.raw_vy**2,axis=0)
            self.sqr_vz = np.sum(self.raw_vz**2,axis=0)
            self.raw_v2 = self.raw_vx**2+self.raw_vy**2+self.raw_vz**2
            self.rel_vx = self.raw_vx-self.mean_vx
            self.rel_vy = self.raw_vy-self.mean_vy
            self.rel_vz = self.raw_vz-self.mean_vz
            self.rel_v2 = self.rel_vx**2+self.rel_vy**2+self.rel_vz**2
            self.cov_v2 = (self.raw_vx-self.mean_vx)**2+\
                          (self.raw_vy-self.mean_vy)**2+\
                          (self.raw_vx-self.mean_vz)**2
            self.rx_hat = np.zeros_like(self.rx_rel)
            self.ry_hat = np.zeros_like(self.rx_rel)
            self.rz_hat = np.zeros_like(self.rx_rel)
            ok = self.r > 0
            self.rx_hat[ok]= self.rx_rel[ok]/self.r[ok]
            self.ry_hat[ok]= self.ry_rel[ok]/self.r[ok]
            self.rz_hat[ok]= self.rz_rel[ok]/self.r[ok]
            self.vr_raw = self.rx_hat*self.raw_vx+\
                      self.ry_hat*self.raw_vy+\
                      self.rz_hat*self.raw_vz
            self.vr_rel = self.rx_hat*self.rel_vx+\
                      self.ry_hat*self.rel_vy+\
                      self.rz_hat*self.rel_vz
            self.vt2_raw = (self.raw_vx-self.vr_raw*self.rx_hat)**2+\
                       (self.raw_vy-self.vr_raw*self.ry_hat)**2+\
                       (self.raw_vz-self.vr_raw*self.rz_hat)**2
            self.vt2_rel = (self.rel_vx-self.vr_rel*self.rx_hat)**2+\
                       (self.rel_vy-self.vr_rel*self.ry_hat)**2+\
                       (self.rel_vz-self.vr_rel*self.rz_hat)**2
            self.zero_raw = self.vt2_raw + self.vr_raw**2 - self.raw_v2
            self.zero_rel = self.vt2_rel + self.vr_rel**2 - self.rel_v2


        self.axis = axis
        if self.axis == 0:
            self.this_h = self.this_y
            self.h_label='y'
            self.this_v = self.this_z
            self.v_label='z'
        elif self.axis == 1:
            self.this_h = self.this_z
            self.h_label='z'
            self.this_v = self.this_x
            self.v_label='x'
        if self.axis == 2:
            self.this_h = self.this_x
            self.h_label='x'
            self.this_v = self.this_y
            self.v_label='y'
