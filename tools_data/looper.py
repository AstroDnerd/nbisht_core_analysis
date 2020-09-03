
import matplotlib
matplotlib.use('Agg')
import yt
import math
import matplotlib.pyplot as plt
import numpy as np
nar = np.array
import time
from importlib import reload
import pyximport; pyximport.install()
import particle_ops
import particle_grid_mask
import loop_tools
reload(loop_tools)
import h5py
#from yt.analysis_modules.level_sets.api import *
from yt.data_objects.level_sets import *
import copy 
import pdb
from importlib import reload
from collections import defaultdict
import weakref
import trackage
import tracks_read_write
import os

def count_particles(fname='tools_data/n_particles.txt'):
    fptr = open(fname,'r')
    lines=fptr.readlines()
    fptr.close()
    parts = np.zeros([len(lines),2])
    pd = {}
    for n,line in enumerate(lines):
        parts[n] = np.array(line.split(),dtype='int')
    pd = dict( zip( parts[:,0], parts[:,1]))
    return pd

def get_all_nonzero(fname='tools_data/n_particles.txt'):
    fptr = open(fname,'r')
    lines=fptr.readlines()
    fptr.close()
    parts = np.zeros([len(lines),2])
    for n,line in enumerate(lines):
        parts[n] = np.array(line.split(),dtype='int')
    all_nonzero = parts[:,0][ parts[:,1] >0]
    core_list = all_nonzero.astype('int')[::-1]
    return core_list #parts, all_nonzero, 
            
class core_looper():
    """We want to analyze the core pre-images for a bunch of frames.
    *core_looper* keeps track of stuff (where the data is, how to loop, all the stuff)
    *snapshot* is for one core in one frame, and packs up the particle values.

    Usage:
        Make a core_looper with a 
         data_directory (where the stuff is)
         data_template  (how to get data from frames)
         out_prefix     (what to call your stuff)
         frame_list     (what frames to analyze)
         target_frame   (what frame the end states are in)
         
        and some extra stuff,
         fields_from_grid (in addition to density and cell volume)
         individual_particle_tracks (doesn't work),

        for example
         >>> u05_loop = looper.core_looper(directory='my_scratch',out_prefix='u05', frame_list=[0,10])

        You need to get the target indices with get_target_indices, this probably needs
        to be worked out better; right now pulling the peaks from an hdf5 file is best.
         >>> u05_loop.get_target_indices(h5name = 'u05_0125_peaklist.h5', core_list=[10])

        Then make some analysis loops, like the ones in run_loop.py
         
    """
    def __init__(self,savefile=None,
                 directory="./", data_template = "%s/DD%04d/data%04d", sim_name='sim',
                 out_prefix="",
                 frame_list=[], core_list=[], target_frame=0,
                 fields_from_grid = [], 
                 individual_particle_tracks=False,
                 derived=[]):
        #set defaults and/or arguments.
        self.current_frame = None
        self.data_template = data_template
        self.sim_name      = sim_name
        self.frame_list     = frame_list
        self.core_list     = core_list
        self.directory     = directory
        self.target_frame  = target_frame
        self.out_prefix    = out_prefix
        self.fields_from_grid = ['density', 'cell_volume'] + fields_from_grid

        #this is not used.
        self.individual_particle_tracks=individual_particle_tracks

        #the track manager.
        self.tr = None

        #defaults for things to be set later
        self.target_indices = {}
        self.ds = None
        self.field_values=None
        self.snaps = defaultdict(dict) 
        #    defaultdict(whatev) is a dict, but makes a new (whatev) by default
        self.ds_list={}
        self.all_data={}
        self.derived=derived

        self.shift = True

        if savefile is not None:
            if not os.path.exists(savefile):
                print("No such file "+savefile)
            else:
                tracks_read_write.load_loop(self,savefile)
                self.tr.sort_time()

        #read from save file.
        #if savefile is not None:
        #    fptr = open(savefile,'r')
        #    lines = fptr.readlines()
        #    fptr.close()
        #    for line in lines:
        #        line = line.strip() #clean off white space
        #        #treat each line as a command to run on 'self'
        #        exec("self.%s"%sline)

    def save(self,fname = "TEST.h5"):
        tracks_read_write.save_loop(self,fname)
    def load_loop(self,fname = "TEST.h5"):
        tracks_read_write.load_loop(self,fname)

    def get_current_frame(self):
        if self.current_frame is None:
            self.current_frame = self.frame_list[0]
        return self.current_frame

    def load(self,frame=None,dummy=False,derived=None):
        """runs yt.load, and saves the ds so we don't have multiples for many cores."""
        if dummy:
            self.ds = None
            self.ds_list[frame]=None
            return None
        if frame is None:
            frame = self.get_current_frame()
        self.filename = self.data_template%(self.directory,frame,frame)
        new_ds = True
        if frame in self.ds_list:
            if self.ds_list[frame] is not None:
                self.ds = self.ds_list[frame]
                new_ds = False
        if new_ds:
            self.ds = yt.load(self.filename)
            if derived is None:
                derived = self.derived
            for add_derived_field in derived:
                add_derived_field(self.ds)
            self.ds_list[frame] = self.ds
        if True:
            self.ds.periodicity = (True,True,True)
        return self.ds

    def get_target_indices(self,target_frame=None,core_list=None,h5_name=None, peak_radius=1.5,
                          bad_particle_list=None):
        if target_frame is None:
            target_frame = self.target_frame
        if core_list is None and self.core_list is not None:
            core_list = self.core_list
        target_ds = self.load(target_frame)
        print("Stuff: ds %s target frame %s cores %s"%(str(target_ds), str(target_frame), str(core_list)))
        new_indices = loop_tools.get_leaf_indices(target_ds,h5_name = h5_name, 
                                     subset = core_list, peak_radius=peak_radius,
                                                 bad_particle_list=bad_particle_list)
        for core_id in new_indices:
            self.target_indices[core_id] = new_indices[core_id]
        

    def get_region(self,frame=None):
        if frame is not None or self.ds is None:
            ds = self.load(frame)
        region = self.ds.all_data()
        return region

    def make_snapshot(self,frame,core_id,dummy_ds=False):
        if core_id in self.snaps[frame]:
            this_snap = self.snaps[frame][core_id]
        else:
            this_snap = snapshot(self,frame,core_id,dummy_ds=dummy_ds)
            self.snaps[frame][core_id] = this_snap # not a weak ref, needs to persist.weakref.proxy(this_snap)
        return this_snap

    def get_tracks(self):
        if self.tr is None:
            self.tr = trackage.track_manager(self)
        for frame in self.frame_list:
            for core_id in self.core_list:
                this_snapshot = self.make_snapshot(frame,core_id)
                if this_snapshot.R_centroid is None:
                    this_snapshot.get_all_properties()
                this_snapshot.get_particle_values_from_grid()
                self.tr.ingest(this_snapshot)
"""
        #this is not used.
        self.individual_particle_tracks=individual_particle_tracks

        #the track manager.
        self.tr = None

        #defaults for things to be set later
        self.target_indices = {}
        self.ds = None
        self.field_values=None
        self.snaps = defaultdict(dict) 
        #    defaultdict(whatev) is a dict, but makes a new (whatev) by default
        self.ds_list={}
"""



class snapshot():
    """For one core and one time, collect the particle positions and whatnot.
    """
    def __init__(self,loop,frame,core_id,dummy_ds=False):
        self.loop           = weakref.proxy(loop) #cyclic references are bad, weakref helps.
        self.target_indices = weakref.proxy(loop.target_indices[core_id])
        self.dummy_ds=dummy_ds
        if dummy_ds:
            self.ds=None
            self.time = -1
        else:
            self.ds             = weakref.proxy(loop.load(frame,dummy=dummy_ds) )
            self.time = self.ds['InitialTime']
        self.core_id        = core_id
        self.frame          = frame

        #stubs for the quantities we'll compute.
        self.mask=None
        self.pos=None
        self.R_centroid=None
        self.field_values={}
        self.R_centroid  =None #(centroid weighted by grid quantities)
        self.R_vec       =None #(particle position relative to centroid)
        self.R_mag       =None #(magnitude of position)
        self.N_vec       =None #(normal vector)
        self.V_bulk      =None #(mean motion)
        self.V_rel       =None #(relative motion)
        self.V_radial    =None #(radial coordinate of velocity)
        
        self.dummy_ds = dummy_ds

    def get_ds(self,frame=None):
        if frame is None:
            frame = self.frame
        self.ds             = weakref.proxy(self.loop.load(frame,dummy=False) )
        return self.ds

    def get_all_properties(self):
        """Run all the relevant analysis pieces."""
        print("get data core %d frame %d"%(self.core_id,self.frame))
        self.get_current_mask()
        self.get_current_pos_vel()
        self.get_particle_values_from_grid()
        self.compute_relative_coords()

    def get_region(self,region_stuff=None):
        """I will probably extend this to make more complex regions."""
        region = self.ds.all_data()
        return region

    def get_current_mask(self):
        """get the particle mask that relates particles for this core_id and this frame
        to the particles in the target_indices from the target_frame"""
        these_pids=self.target_indices.astype('int64')
        if type(these_pids) == yt.units.yt_array:
            these_pids = these_pids.v
        data_region = self.get_region(self.frame)
        mask_to_get=np.zeros(these_pids.shape,dtype='int32')
        my_indices = data_region['particle_index'].astype('int64')
        found_any, mask = particle_ops.mask_particles(these_pids,my_indices,mask_to_get)
        self.mask = mask
        return found_any, mask
        
    def check_and_fix_bad_particles(self,mask):
        if mask.sum() == self.target_indices.size:
            return
        else:
            print("WARNING:  missing particle. Likely out of grid bounds.")
        data_region = self.get_region(self.frame)
        these_pids=self.target_indices.astype('int64')
        my_indices = data_region['particle_index'].astype('int64')
        bad_particles = []
        for ppp in these_pids:
            if ppp not in my_indices:
                bad_particles.append(ppp)
        if len(bad_particles) == 0:
            print("WORSE WARNING: still can't find the missing particles.")
            return
        grids = []
        found_particles=[]
        for grid in self.ds.index.grids:
            for ppp in bad_particles:
                if ppp in grid['particle_index']:
                    grids.append(grid)
                    found_particles.append(ppp)
                    ok = np.where( grid['particle_index'] == ppp)[0]
                    self.pos = self.ds.arr( np.concatenate([self.pos,grid['particle_position'][ok]]),'code_length')
                    self.vel = self.ds.arr( np.concatenate([self.vel,grid['particle_velocity'][ok]]),'code_velocity')
                    self.ind = self.ds.arr( np.concatenate([self.ind,grid['particle_index'][ok]]),'dimensionless')

        if len(found_particles) == 0:
            print("EVEN WORSE WARNING: still can't find the missing particles.")
            return
        if len(found_particles) != len(bad_particles):
            print("WARNING: found some but not all of the bad particles.")







    def get_current_pos_vel(self):

        if self.mask is not None:
            found_any, mask = self.get_current_mask()
            if not found_any:
                raise
        region = self.get_region(self.frame)
        self.ind = region['particle_index'][mask == 1]
        args=np.argsort(self.ind)
        self.ind=self.ind[args]
        self.pos = region['particle_position'][mask == 1][args]
        self.vel = region['particle_velocity'][mask == 1][args]
        #self.check_and_fix_bad_particles(mask)


    def compute_relative_coords(self):
        """Compute and store 
        R_centroid (centroid weighted by grid quantities)
        R_vec      (particle position relative to centroid)
        R_mag      (magnitude of position)
        N_vec      (normal vector)
        V_bulk     (mean motion)
        V_rel      (relative motion)
        V_radial   (radial coordinate of velocity)
        Assumes that the following have been set properly:
        self.field_values
        self.pos
        self.vel
        self.ds
        """
        if len(self.field_values.keys())==0:
            self.get_particle_values_from_grid()

        m = self.field_values['density'].sum()
        if self.loop.shift:
            shifted = loop_tools.shift_particles(self.ds,self.pos,shiftRight=False)
        else:
            shifted = copy.copy(self.pos)
        if self.ds is not None:
            self.R_centroid = self.ds.arr([0,0,0],'code_length')
            self.V_bulk = self.ds.arr([0,0,0],'code_velocity')
        else:
            self.R_centroid = yt.units.yt_array([0,0,0],'cm')
            self.V_bulk =    yt.units.yt_array([0,0,0],'cm/s')
        for dim in range(3):
            self.R_centroid[dim] = (shifted[:,dim]*self.field_values['density']).sum()/m
            self.V_bulk[dim] = (self.vel[:,dim]*self.field_values['density']).sum()/m
        self.R_vec = shifted - self.R_centroid
        self.R_mag = (self.R_vec**2).sum(axis=1)**0.5
        self.N_vec = np.zeros_like(self.R_vec)
        for dim in range(3):
            self.N_vec[:,dim] = self.R_vec[:,dim]/self.R_mag
        self.V_relative = self.vel - self.V_bulk
        self.V_radial = (self.V_relative * self.N_vec).sum(axis=1)
        self.field_values['V_radial']=self.V_radial
        #self.field_values['V_mag']=((self.V_relative*self.V_relative).sum())**0.5
        self.field_values['R_mag']=self.R_mag
        #self.field_values['V_radial']=self.V_radial
        #self.field_values['V_relative']=self.V_relative


    def get_particle_values_from_grid(self, field_list=[]):
        """Get the simple nearest-sampled grid point from the particles in self.pos.
        Assumes that self.pos has been set correctly.
        """
        if self.pos is None:
            self.get_current_pos_vel()

        if self.field_values is None:
            self.field_values={}
        fields_to_get = np.unique(self.loop.fields_from_grid + field_list)
        number_of_new_fields=0

        for field in fields_to_get:
            #initialize field values to -1.
            if field not in self.field_values:
                self.field_values[field] = np.zeros(self.ind.size)-1
                number_of_new_fields += 1

        if number_of_new_fields == 0:
            return

        good_index_sort_np = np.array(copy.copy(self.ind)).astype('int64')
        for grid in self.ds.index.grids[-1::-1]:
            grid_selector = np.zeros([3,good_index_sort_np.size],dtype='int32') #grid i,j,k index selector.
            particle_selector = np.zeros(good_index_sort_np.size,dtype='int32') #mask between all particles and this grid.
            mask_to_get_3 = np.zeros(good_index_sort_np.shape, dtype='int32')   #particles in this grid.  This is only for existence checking.
            #this gives the particles that live in _this grid._
            found_any_g, mask_g = particle_ops.mask_particles(good_index_sort_np, grid['particle_index'].astype('int64'), mask_to_get_3)
            if found_any_g:
                particle_grid_mask.particle_grid_mask_go_i2(self.pos[:,0],self.pos[:,1],self.pos[:,2], 
                                                            grid.LeftEdge, grid.dds, 
                                                            grid.ActiveDimensions,grid.child_mask, 
                                                            grid_selector,
                                                            particle_selector)

                particle_selector = particle_selector==1
                if particle_selector.sum() == 0:
                    continue
                for field in self.field_values:
                    values = grid[field][[grid_selector[i][particle_selector] for i in [0,1,2]]]
                    self.field_values[field][particle_selector] = values
                if self.loop.individual_particle_tracks:
                    print("ERROR: individual particle tracks not implimented")
                    raise
                    particles_in_grids = these_pids[particle_selector]
                    for p in particles_in_grids:
                       pdict[p]=pdict.get(p,{})
                       for k in ['times','cycles']:
                           if not pdict[p].has_key(k): pdict[p][k]=[]
                       pdict[p]['times'].append(times[-1])
                       pdict[p]['cycles'].append(cycles[-1])
                    for field in ['density']:
                       values = grid[field][[grid_selector[i][particle_selector] for i in [0,1,2]]]
                       for n,p in enumerate(particles_in_grids):
                           pdict[p][field] = pdict[p].get(field,[])
                           pdict[p][field].append(values[n])
                          
        #error check.
        if  (self.field_values['cell_volume'] < 0).any():
            NM= (self.field_values['cell_volume'] < 0).sum()
            print("ERROR: some particles (%d of them) not found.  This is problematic."%NM)

#Decorators frame_loop and core_loop takes functions and makes them loop over frame and core_id.
#The decorator takes care of loading the data, so your funciton only needs to use it.
#
#1.) write a function that does stuff as a method of the clump_looper
#2.) when you define it, put @looper_loop above the definition.
#3.) add it to the function
#4.) looper_loop shouldn't be modified
#See the 'full_plot' definition below,
#Also See decorator_test.py for more details.
def frame_loop(function):
    def wrapper(self,*args,**kwargs):
        for self.current_frame in self.frame_list:
            self.ds = self.load(self.current_frame)
            function(self,*args,**kwargs)
    return wrapper

def core_loop(function):
    def wrapper(looper,*args,**kwargs):
        for looper.current_frame in looper.frame_list:
            snapshot.ds = weakref.proxy(looper.load(frame=looper.current_frame))
            for looper.core_id in looper.core_list:
                this_snapshot = looper.make_snapshot(looper.current_frame,looper.core_id)
                if this_snapshot.R_centroid is None:
                    this_snapshot.get_all_properties()
                output = function(looper,this_snapshot,*args,**kwargs)
        return output #this is what the function will return
    return wrapper    #this return is from the decorator
# Decorators?
# A decorator looks like
# @some_decorator
# def some_function(arguments):
#    do_stuff()
# The decorator then does stuff to the way some_function behaves.
# So, if my decorator is do_five_times,
# def do_five_times(function):
#    def wrapper():
#        for i in range(5):
#            function()
# then I do
# @do_five_times
# def say_no():
#    print("no")
#
# say_no()
# will get called five times.  
# The decorator does
# say_no = do_five_times(say_no)
