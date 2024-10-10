from starter2 import *
reload(track_info)


#
# Define tracks, including track location and build information.
# The constructor also registers each one in track_info.tracks
#

track_info.track('nb101',
                 sim_directory='/data/cb1/nbisht/anvil_scratch/projects/128/B2',
                 target_frame = 125,
                 mountain_top = "/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_mountain_tops.h5",
                 peak_fname = '/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_0125_peaklist.h5',
                 track_file = "/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_all_frames.h5",
                 bad_particles="/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_bad_particles.h5",
                 mode_fname = "/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/core_formation_mode_new_nb101.h5",
                 export_to_ML_fname = "/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_ML_dataset.csv",
                 SinkClumpLink_fname = "/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/NonsinkSimClump_SinkSimSinks_link.json",
                 frame_list = "all_frames",
                 field_list = track_info.field_lists['primitive'])

track_info.track('nb102',
                 sim_directory='/data/cb1/nbisht/anvil_scratch/projects/128/B2_sink_nazare',
                 target_frame = 125,
                 mountain_top = "/data/cb1/nbisht/anvil_scratch/projects/128/B2_sink_nazare/datasets/nb102_mountain_tops.h5",
                 peak_fname = '/data/cb1/nbisht/anvil_scratch/projects/128/B2_sink_nazare/datasets/nb102_0125_peaklist.h5',
                 track_file = "/data/cb1/nbisht/anvil_scratch/projects/128/B2_sink_nazare/datasets/nb102_every_10.h5",
                 bad_particles="/data/cb1/nbisht/anvil_scratch/projects/128/B2_sink_nazare/datasets/nb102_bad_particles.h5",
                 mode_fname = "/data/cb1/nbisht/anvil_scratch/projects/128/B2_sink_nazare/datasets/core_formation_mode_new_nb102.h5",
                 export_to_ML_fname = "/data/cb1/nbisht/anvil_scratch/projects/128/B2_sink_nazare/datasets/nb102_ML_dataset.csv",
                 frame_list = "every_ten",
                 field_list = track_info.field_lists['primitive'])


