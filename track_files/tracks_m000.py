
from starter2 import *
reload(track_info)

# m_ for multiple target frames
# track name: m02targetframe, where 02 is the second simulation of the 200 series

#
# Define tracks, including track location and build information.
# The constructor also registers each one in track_info.tracks
#

track_info.track('m0250',
                 sim_directory=dl.sims['u202'],
                 target_frame = 50, 
                 mountain_top = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/SFR_game/m0250_mountain_top_50.h5",
                 peak_fname = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/SFR_game/m0250_peaklist_50.h5",
                 track_file = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/SFR_game/m0250_every10_50.h5",
                 mode_fname = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/SFR_game/m0250_coreformationmode_50.h5",
                 frame_list = "every_ten",  
                 field_list = track_info.field_lists['minimal'])

track_info.track('m0260',
                 sim_directory=dl.sims['u202'],
                 target_frame = 60,
                 mountain_top = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/SFR_game/m0260_mountain_top_60.h5",
                 peak_fname = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/SFR_game/m0260_peaklist_60.h5",
                 track_file = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/SFR_game/m0260_every10_60.h5",
                 mode_fname = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/SFR_game/m0260_coreformationmode_60.h5",
                 frame_list = "every_ten",  
                 field_list = track_info.field_lists['minimal'])

track_info.track('m0270',
                 sim_directory=dl.sims['u202'],
                 target_frame = 70,
                 mountain_top = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/SFR_game/m0270_mountain_top_70.h5",
                 peak_fname = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/SFR_game/m0270_peaklist_70.h5",
                 track_file = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/SFR_game/m0270_every10_70.h5",
                 mode_fname = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/SFR_game/m0270_coreformationmode_70.h5",
                 frame_list = "every_ten",  
                 field_list = track_info.field_lists['minimal'])

track_info.track('m0280',
                 sim_directory=dl.sims['u202'],
                 target_frame = 80,
                 mountain_top = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/SFR_game/m0280_mountain_top_80.h5",
                 peak_fname = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/SFR_game/m0280_peaklist_80.h5",
                 track_file = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/SFR_game/m0280_every10_80.h5",
                 mode_fname = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/SFR_game/m0280_coreformationmode_j0.h5",
                 frame_list = "every_ten",  
                 field_list = track_info.field_lists['minimal'])
