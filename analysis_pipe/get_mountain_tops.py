from starter2 import *

import mountain_top
reload(mountain_top)

import coreset_data
reload(coreset_data)

#
# Mountain top removal tool.
# 1.) Find all peaks with rho_max > 1e4 (set in coreset_data with the verifier)
# 2.) Check for overlap between different cores.
# 3.) Separate peaks by contouring
# 4.) Re-run mountain top removal.

def get_mountain_tops(trackname):
    this_track = track_info.tracks[trackname]
    mountain_top_name = this_track.mountain_top

    if os.path.exists(mountain_top_name):
        print("File exists, not performing mountain tops.", mountain_top_name)
        return 0

    do_mountain_projections=True

    verifier = coreset_data.verify_cores_generic
    this_radius_dict={}

#
# 1.) Make all the mountain tops.  First pass.
#
    leaf_storage_tmp={}
    kludge={}
    #kludge={}
    MT=mountain_top.cut_mountain_top( trackname, target_fname = mountain_top_name, 
                                     do_projections=False, verify=verifier,
                                     leaf_storage=leaf_storage_tmp, kludge=kludge, 
                                     radius_dict=this_radius_dict)
    leaf_storage_all = leaf_storage_tmp


#
# 2.) Check for overlap.
#
    overlap_all=mountain_top.check_overlap(leaf_storage_all)


#
# 3.) Separate cores with overlap.
#

    new_thresholds = {}  #this gets filled with new thresholds
    if 'all_peaks_tmp' not in dir():
        all_peaks_tmp=defaultdict(dict)     #this is just to reduce redundant work when debugging.
    mountain_top.split_all(trackname, overlap_all, new_thresholds, all_peaks_tmp, 
                           do_projections=False, radius_dict=this_radius_dict)

#
# 4.) Make all the leaves again, with new thresholds.  Makes PEAK images of all peaks, even rejects.
#
    leaf_storage_2={}
    kludge={}
    MT=mountain_top.cut_mountain_top( trackname, target_fname = mountain_top_name, 
                                     do_projections=True, verify=verifier,
                                     kludge=kludge, leaf_storage=leaf_storage_2, 
                                     cut_override = new_thresholds, 
                                     radius_dict=this_radius_dict)
    overlap_2=mountain_top.check_overlap(leaf_storage_2)
