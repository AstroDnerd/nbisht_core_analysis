
from starter2 import *
import davetools
reload(davetools)
import p49_fields
reload(p49_fields)
import math
import pcolormesh_helper as pch
from matplotlib.ticker import PercentFormatter
# --- --- --- --- --- --- ---

class visit(): 
    def __init__(self,the_loop):
        self.this_looper = the_loop
        self.cores_used = []
        self.bad_cores = []

    def qtyRun(self,sim,core_list=None): 
        print('inside qtyRun')
        thtr = self.this_looper.tr

        the_level = 4
        # CORES
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        # THE FRAMES
        the_frame = thtr.frames[-1:] 

        # CORE-LOOP
        for nc,core_id in enumerate(core_list):
            self.cores_used.append(core_id)

            ds = self.this_looper.load(the_frame[0]) 
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            ms.particle_pos(core_id)

            the_center = ms.mean_center[:,-1]
            the_radius = 1/128

            # MAKING SPHERES
            if 0:
                the_sphere = ds.sphere(the_center, the_radius)
                rho = the_sphere['density']
                bx = the_sphere['magnetic_field_x']
                by = the_sphere['magnetic_field_y']
                bz = the_sphere['magnetic_field_z']
                vx = the_sphere['velocity_x']
                vy = the_sphere['velocity_y']
                vz = the_sphere['velocity_z']
                
            # MAKING COVERING GRIDS
            if core_id == 228:
                print(core_id)
                delta_r = the_center + (the_radius/2)
                delta_l = the_center - (the_radius/2)
                the_redge = ds.arr(delta_r,'code_length')
                the_ledge = ds.arr(delta_l,'code_length')
                min_dx = ds.index.get_smallest_dx()
                the_size = ((the_redge - the_ledge)/min_dx).astype('int')
                cg = ds.covering_grid(the_level,the_ledge,the_size)
               
                rho = cg['density'] 
                bx = cg['magnetic_field_x']
                by = cg['magnetic_field_y']
                bz = cg['magnetic_field_z']
                vm = cg['velocity_magnitude']
                vx = cg['velocity_x']
                vy = cg['velocity_y']
                vz = cg['velocity_z']


        # CACHE TO DISK 
        if 1: 
            hfivename = 'p66_brho/h5files/covergrid_%s_228_two.h5'%(sim)#,core_id)
            Fptr = h5py.File(hfivename,'w')
            
            Fptr['bx']=bx
            Fptr['by']=by
            Fptr['bz']=bz
            Fptr['rho']=rho
            Fptr['velocity_mag']=vm
            Fptr['velocity_x']=vx
            Fptr['velocity_y']=vy
            Fptr['velocity_z']=vz


# MAIN
import three_loopers_six as TL6
if 'clobber' not in dir():
    clobber=True
if 'scope1' not in dir() or clobber:
    scope1=visit(TL6.loops['u601'])
    core_list1 = TL6.loops['u601'].core_by_mode['Alone']
if 'scope2' not in dir() or clobber:
    scope2=visit(TL6.loops['u602'])
    core_list2 = TL6.loops['u602'].core_by_mode['Alone']
if 'scope3' not in dir() or clobber:
    scope3=visit(TL6.loops['u603'])
    core_list3 = TL6.loops['u603'].core_by_mode['Alone']

simnames = ['u601','u602', 'u603']
for nt,tool in enumerate([scope1]):#,scope2,scope3]):
    # WHICH CORES & RUN IF NOT CACHED
    if 1:
        all_cores = np.unique(tool.this_looper.tr.core_ids)
        #core_list = all_cores[:1]  #DEBUG
        core_list = all_cores
        #core_list = core_list1 
        tool.qtyRun(nt,core_list=core_list) 



