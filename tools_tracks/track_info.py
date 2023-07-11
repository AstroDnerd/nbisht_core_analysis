from starter2 import *
import xtra_energy

#
# The track information object.  
# Contains info about a track, including build information.
# Also contains typical field lists, for consistency.
#

if 'tracks' not in dir():
    tracks={}

field_lists={}
field_lists['minimal']=  [YT_x,YT_y,YT_z,YT_density, YT_cell_volume]
field_lists['primitive']=field_lists['minimal']+\
        [YT_magnetic_field_x, YT_magnetic_field_y, YT_magnetic_field_z]+\
        [YT_velocity_x,YT_velocity_y,YT_velocity_z] +\
        [YT_velocity_magnitude,YT_magnetic_field_strength, YT_velocity_divergence] +\
        [YT_potential_field]
field_lists['acceleration'] =  [YT_acceleration_x, YT_acceleration_y, YT_acceleration_z]
field_lists['particle_pos'] =  [YT_particle_pos_x, YT_particle_pos_y, YT_particle_pos_z, YT_particle_index]
field_lists['shear_tensor']= [YT_dxvx,YT_dxvy,YT_dxvz]+\
                             [YT_dyvx,YT_dyvy,YT_dyvz]+\
                             [YT_dzvx,YT_dzvy,YT_dzvz]
field_lists['most_things'] = field_lists['primitive']+field_lists['acceleration']+field_lists['particle_pos']
field_lists['everything'] = field_lists['most_things']+field_lists['shear_tensor']

#
# Standard Derived field list
#
derived=[xtra_energy.add_force_terms, xtra_energy.add_v_grad]


class track():
    def __init__(self,name,sim_directory=None,target_frame=None,mountain_top=None,peak_fname=None,track_file=None, field_list=None, frame_list=None, clump_parameters=None,
                bad_particles=None, plot_directory="./plots_to_sort", derived_fields=None, mode_fname=None):
        self.name = name
        self.sim_directory =sim_directory
        self.target_frame  =target_frame
        self.mountain_top  =mountain_top
        self.peak_fname     =peak_fname
        self.track_file    =track_file
        self.field_list    = field_list
        self.mode_fname    = mode_fname
        if frame_list is not None:
            if type(frame_list)==str:
                self.frame_list = self.make_frame_list(frame_list)
            else:
                self.frame_list=frame_list
        if clump_parameters is not None:
            self.clump_parameters = clump_parameters
        else:
            #Python does funny things if you make default parameters as mutable objects.
            self.clump_parameters = {}
        self.bad_particles=bad_particles
        tracks[name]=self
        if derived_fields == None:
            self.derived_fields = derived
        else:
            self.derived_fields = derived_fields
        self.plot_directory=plot_directory
    def make_frame_list(self,frame_list):
        target_frame=self.target_frame
        if frame_list ==  'all_frames':
            frame_list = list(range(0,target_frame+1))
        elif frame_list == 'first_4':
            frame_list = list(range(0,3))
        elif frame_list == 'every_ten':
            frame_list = list(range(0,target_frame,10)) + [target_frame]
        else:
            pdb.set_trace()
            print("Bad frame list in track.")
        return frame_list

