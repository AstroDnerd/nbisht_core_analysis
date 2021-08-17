
from starter2 import *
import xtra_energy
import data_locations as dl
reload(dl)
#
# set sim
#
if 'this_simname' not in dir():
    this_simname = 'u11'

for this_simname in ['u201','u202','u203']:
    if this_simname not in ['u201']:
        continue


    output_name = '%s_first_last_t3_nXXX0.h5'%(this_simname)
    if 'this_looper' not in dir():
        #
        # make first_last object
        #
        core_list =  [0]
        target = dl.target_frames[this_simname]
        frame_list = [0, target]
        frame_list = list(range(0,target,10)) + [target]
        frame_list = [100]

        #frame_list = list(range(0, dl.target_frames[this_simname],10))+[dl.target_frames[this_simname]]
        fields = ['x','y','z','density','velocity_magnitude','magnetic_field_strength', 'velocity_divergence']
        fields += ['velocity_x','velocity_y','velocity_z']
        fields += ['magnetic_field_%s'%s for s in 'xyz']
        fields += ['PotentialField']
        output_base = "primitive_test"
        derived=[]


        this_looper = looper.core_looper(directory= dl.sims[this_simname],
                                         sim_name = this_simname,
                                         out_prefix = this_simname,
                                         target_frame = dl.target_frames[this_simname],
                                         frame_list = frame_list,
                                         core_list =   core_list,
                                         fields_from_grid=fields,
                                         derived = derived
                                      )
        ds = this_looper.load(frame_list[0])
        ad = ds.all_data() #ds.region([0.5]*3,[0.4]*3,[0.6]*3)
        this_looper.target_indices[0]=ad['particle_index']
        this_looper.get_tracks()
        this_looper.save(output_name)
