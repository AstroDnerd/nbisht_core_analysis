from starter2 import *

#First install track info into one of the track modules.
# e.g. tools_data/tracks_t200.py
# Can be a new tracks_whatever.py, but must be imported by track_loader.py
# Then run this.

this_trackname = 't002'
import get_peaks
get_peaks.get_peaks(this_trackname)
import image_peaks
#image_peaks.image_peaks(this_trackname)
import get_mountain_tops
get_mountain_tops.get_mountain_tops(this_trackname)
import new_tracks
new_tracks.get_tracks(this_trackname)
import image_mountain
#image_mountain.image_mountains(this_trackname)
import mode_grader
reload(mode_grader)
delta_alone = 0.025
#mode_grader.grade_modes(this_trackname, delta_alone=delta_alone)
import buddy_centroid
reload(buddy_centroid)
#buddy_centroid.buddy_centroid(this_trackname, thresh=2*delta_alone)
import rho_time
#rho_time.run(this_trackname)
import browser_skin
reload(browser_skin)
browser_skin.make_browser(this_trackname)
