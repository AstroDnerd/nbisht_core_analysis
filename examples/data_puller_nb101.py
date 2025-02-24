"""
This just pulls partilce data from enzo and stores it.
Changing
core_list 
frame_list
fields
changes what gets extracted.
"""
from starter2 import *
import xtra_energy
import data_locations as dl
reload(dl)
import looper
#
# set sim
#
if 'this_simname' not in dir():
    this_simname = 'nb101'
all_nonzero = looper.get_all_nonzero(dl.n_particles[this_simname])

output_base = "%s_cores"%this_simname
if 0:
    """this set of parameters extracts all primitive quantities"""
    core_list =  all_nonzero
    frame_list = [0]#[0]#+list(range(10,130,10))+[125]
    fields = ['x','y','z','density','velocity_magnitude','magnetic_field_strength', 'velocity_divergence']
    fields += ['velocity_x','velocity_y','velocity_z']
    fields += ['magnetic_field_%s'%s for s in 'xyz']
    fields += ['PotentialField']
    output_base = "primitive_test"
    derived=[]

if 0:
    """this set of parameters extracts all primitive quantities"""
    core_list = all_nonzero.astype('int')[::-1][3:4]
    target_frame = dl.target_frames[this_simname]
    frame_list = [0]# [0]+list(range(10,target_frame,10))+[target_frame]
    fields = ['x','y','z','density']
    output_base = "%s_density_only"%this_simname
    derived=[]

if 0:
    """this set of parameters extracts all primitive quantities"""
    core_list = all_nonzero.astype('int')[::-1]
    target_frame = dl.target_frames[this_simname]
    frame_list = list(range(0,target_frame,10))+[target_frame]
    fields = ['x','y','z','velocity_magnitude','magnetic_field_strength', 'velocity_divergence']
    fields += ['velocity_x','velocity_y','velocity_z']
    fields += ['magnetic_field_%s'%s for s in 'xyz']
    fields += ['PotentialField']
    output_base = "%s_all_primitives"%this_simname
    derived=[]

if 0:
    """This set extracts magnetic work"""
    core_list = all_nonzero.astype('int')
    frame_list =[0,1]+list(range(10,130,50))+[125]
    fields = ['mag_work']
    derived=[xtra_energy.add_force_terms]
    output_base = 'mag_work_only'

if 1:
    """This set extracts magnetic work"""
    core_list = all_nonzero.astype('int')
    target_frame = dl.target_frames[this_simname]
    frame_list =list(range(0,target_frame,10))+[target_frame]
    fields = ['x','y','z','density'] #always need these
    fields += ['kinetic_energy_density','therm_energy','magnetic_energy','PotentialField','vorticity_magnitude']
    fields += ['velocity_x','velocity_y','velocity_z']
    fields += ['magnetic_field_%s'%s for s in 'xyz']
    fields += ['PotentialField']
    derived=[xtra_energy.add_force_terms, xtra_energy.add_energies]
    output_base = '%s_energy_vorticity'%this_simname
    #pdb.set_trace()

if 1:
    #Pull a whole list for each core, saving each core to its own file
    for core in core_list:  
        output_name = '%s_primitives_c%04d_nXXX0.h5'%(output_base,core)
        if os.path.exists(output_name):
            print("File exists, skipping "+output_name)
            continue
        this_looper = looper.core_looper(directory= dl.sims[this_simname],
                                         sim_name = this_simname,
                                         out_prefix = this_simname,
                                         target_frame = dl.target_frames[this_simname],
                                         frame_list = frame_list,
                                         core_list =  [core],# core_list,
                                         fields_from_grid=fields,
                                         derived = None#derived
                                      )
        this_looper.get_target_indices(h5_name=dl.peak_list[this_simname],)
                                         #bad_particle_list=dl.bad_particles.get(this_simname,None))
        this_looper.get_tracks()
        this_looper.save(output_name)

if 0:
    #Pull all cores at once. 
    output_name = '%s_many.h5'%(output_base)
    if not os.path.exists(output_name):
        this_looper = looper.core_looper(directory= dl.sims[this_simname],
                                         sim_name = this_simname,
                                         out_prefix = this_simname,
                                         target_frame = dl.target_frames[this_simname],
                                         frame_list = frame_list,
                                         core_list =  core_list,
                                         fields_from_grid=fields,
                                         derived = derived
                                      )
        this_looper.get_target_indices(h5_name=dl.peak_list[this_simname],
                                         bad_particle_list=dl.bad_particles[this_simname])
        this_looper.get_tracks()
        this_looper.save(output_name)
