from starter2 import *
import track_info
reload(track_info)



#
# Define tracks, including track location and build information.
# The constructor also registers each one in track_info.tracks
#

track_info.track('t202',
                 sim_directory = dl.sims['u202'],
                 target_frame = 60, #half a free fall time
                 peak_fname = 'datasets_small/t202_peaks.h5',
                 mountain_top = 'datasets_small/t202_mountain_top.h5',
                 frame_list = 'every_ten',
                 field_list = track_info.field_lists['primitive'],
                 clump_parameters = {'small_test':False}
                )

