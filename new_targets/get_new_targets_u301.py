from starter2 import *

import mountain_top
reload(mountain_top)
reload(looper)
import tracks_read_write
reload(tracks_read_write)
#
# get mountain tops
#

this_simname = 'u301'
mountain_top_name = "%s_mountain_tops_take_8.h5"%this_simname
do_mountain_projections=True

import coreset_data
reload(coreset_data)
verifiers = {'u302':coreset_data.verify_cores_u302, 'u301':coreset_data.verify_cores_u301,
             'u303':coreset_data.verify_cores_u303}

this_radius_dict={'u301':coreset_data.radius_u301,
                  'u302':coreset_data.radius_u302,
                  'u303':coreset_data.radius_u303}[this_simname]

#
# Make all the leaves and their overlaps.
#
if 'leaf_storage_all' not in dir():#'leaf_storage' not in dir() and False:
    leaf_storage_all={}
    k0,k1=11,14
    #kludge={'peak_id':list(range(100,200))}
    kludge={}
    MT=mountain_top.cut_mountain_top( this_simname, target_fname = mountain_top_name, 
                                     do_projections=do_mountain_projections, verify=verifiers[this_simname],
                                     leaf_storage=leaf_storage_all, kludge=kludge, 
                                     radius_dict=this_radius_dict)
    overlap_all=mountain_top.check_overlap(leaf_storage_all)


if 1:
    new_thresholds = {}  #this gets filled with new thresholds
    if 'all_peaks_tmp' not in dir():
        all_peaks_tmp=defaultdict(dict)     #this is just to reduce redundant work when debugging.
    mountain_top.split_all(this_simname, overlap_all, new_thresholds, all_peaks_tmp, 
                           do_projections=True, radius_dict=this_radius_dict)

#
# Make all the leaves again, with new thresholds
#
if 1:
    leaf_storage_2={}
    #kludge={'peak_id':[k0,k1]}
    kludge={}
    MT=mountain_top.cut_mountain_top( this_simname, target_fname = mountain_top_name, 
                                     do_projections=do_mountain_projections, verify=verifiers[this_simname],
                                     kludge=kludge, leaf_storage=leaf_storage_2, 
                                     cut_override = new_thresholds, 
                                     radius_dict=this_radius_dict)
    overlap_2=mountain_top.check_overlap(leaf_storage_2)
