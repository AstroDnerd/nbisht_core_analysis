from starter2 import *
import track_loader
#
# This is the basic script for the full analysis pipeline.
#
# To Do:
#    1.) Look at track_files/tracks_t000.py
#    2.) Make a new track file in track_files OR add to an existing one.
#    3.) If you make a new one, import it in the track_loader.py
#    4.) Name your tracks something sensible.
#    5.) Run this code.
#    
#

#First install track info into one of the track modules.
# e.g. tools_data/tracks_t200.py
# Can be a new tracks_whatever.py, but must be imported by track_loader.py
# Then run this.

this_trackname = 'nb101'
sink_trackname = 'nb102'

#Use total density (density+sink density) or just density
USE_TOTAL = False
density_placeholder = 'density'
if USE_TOTAL:
    density_placeholder = 'total_density'
    sink_dict_loc = '/data/cb1/nbisht/anvil_scratch/projects/128/B2_sink_nazare/datasets/sink_clump_position_with_sink_totaldensity.json'
    
else:
    sink_dict_loc = track_info.tracks[this_trackname].SinkClumpLink_fname

import get_peaks
print("GET PEAKS")
get_peaks.get_peaks(this_trackname, density_placeholder)

import image_peaks
#image_peaks.image_peaks(this_trackname)

#FOR LINKING NONSINK CORES TO SINKS
#import SinkClumpLink
#SinkClumpLink.SinkClumpLinker(this_trackname, sink_trackname, do_projections = True)

import get_mountain_tops
print("GET MOUNTAIN")
get_mountain_tops.get_mountain_tops(this_trackname, density_placeholder)

import new_tracks
print("GET TRACKS")
new_tracks.get_tracks(this_trackname, density_placeholder)

#import image_mountain
#image_mountain.image_mountains(this_trackname)

import mode_grader
reload(mode_grader)
delta_alone = 0.025
mode_grader.grade_modes(this_trackname, delta_alone=delta_alone)  

#import buddy_centroid
#reload(buddy_centroid)
#new_core_list = buddy_centroid.buddy_centroid(this_trackname, thresh=2*delta_alone, sink_dict_loc=sink_dict_loc)  

#import rho_time  
#rho_time.run(this_trackname, new_core_list)

#import browser_skin
#reload(browser_skin)
#browser_skin.make_browser(this_trackname)
