from starter2 import *

import mountain_top
reload(mountain_top)
reload(looper)
import tracks_read_write
reload(tracks_read_write)
#
# get mountain tops
#

this_simname = 'u703'
mountain_top_name = "%s_mountain_tops_repeat.h5"%this_simname
do_projections=True

import coreset_data
reload(coreset_data)
verifiers = {'u702':coreset_data.verify_cores_u702, 'u701':coreset_data.verify_cores_u701,
             'u703':coreset_data.verify_cores_u703}

this_radius_dict={'u701':coreset_data.radius_u301,
                  'u702':coreset_data.radius_u302,
                  'u703':coreset_data.radius_u303}[this_simname]

#
# Make all the leaves and their overlaps.
#
if 'leaf_storage_all' not in dir():#'leaf_storage' not in dir() and False:
    leaf_storage_all={}
    #kludge={'peak_id':[10]}
    #kludge={'peak_id':[11, 111, 129, 132, 133, 134, 153, 170, 177, 199]}
    kludge={}
    MT=mountain_top.cut_mountain_top( this_simname, target_fname = mountain_top_name, 
                                     do_projections=False, verify=verifiers[this_simname],
                                     leaf_storage=leaf_storage_all, kludge=kludge, 
                                     radius_dict=this_radius_dict)
    overlap_all=mountain_top.check_overlap(leaf_storage_all)


if 1:
    #this gets filled with new thresholds
    new_thresholds = {}  

    
    #this is just to reduce redundant work when debugging.
    if 'all_peaks_tmp' not in dir():
        all_peaks_tmp=defaultdict(dict)     
        all_peaks_tmp[4]={5:740, 6:740,20:740}
        all_peaks_tmp[5]={6:740}
        all_peaks_tmp[6]={5:740}
        all_peaks_tmp[20]={4:300}
        all_peaks_tmp[10]={11:-1,12:-1}
    mountain_top.split_all(this_simname, overlap_all, new_thresholds, all_peaks_tmp, 
                           do_projections=False, radius_dict=this_radius_dict,
                           contour_step=2,
                          auto_double=False)

#
# Make all the leaves again, with new thresholds
#
if 1:
    leaf_storage_2={}
    #kludge={'peak_id':[k0,k1]}
    kludge={}
    MT=mountain_top.cut_mountain_top( this_simname, target_fname = mountain_top_name, 
                                     do_projections=True, verify=verifiers[this_simname],
                                     kludge=kludge, leaf_storage=leaf_storage_2, 
                                     cut_override = new_thresholds, 
                                     radius_dict=this_radius_dict)
    overlap_2=mountain_top.check_overlap(leaf_storage_2)
    print("Overlap ", overlap_2)
