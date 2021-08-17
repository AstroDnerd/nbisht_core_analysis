
from starter2 import *

import mountain_top
reload(mountain_top)
reload(looper)
import tracks_read_write
reload(tracks_read_write)
#
# get mountain tops
#

this_simname = 'u302'
core_id = 14
mountain_top_name = "u302_mountain_overlap_test_2.h5"
do_mountain_projections=False

def verify_cores_u302(leaf, peak_density, peak_id):
    out = True
    if peak_density < 1100:
        out = False
    skip = [4, 89]
    if peak_id in skip:
        out=False
    return out

if 'leaf_storage_all' not in dir():#'leaf_storage' not in dir() and False:
    kludge={}
    leaf_storage_all={}
    MT=mountain_top.cut_mountain_top( this_simname, target_fname = mountain_top_name, do_projections=do_mountain_projections, verify=verify_cores_u302,kludge=kludge, leaf_storage=leaf_storage_all)#,
                                     #density_cut_fname="u302_contour_mins.h5")
    overlap_all=mountain_top.check_overlap(leaf_storage_all)

if 'leaf_storage' not in dir():#'leaf_storage' not in dir() and False:
    #kludge={'peak_id':[0,1,364, 113,u302_was258]})
    #kludge={'peak_id':[258]}
    #kludge={'peak_id':[core_id]}
    kludge={}
    #kludge={'peak_id':[89,11,14,350, 349, 368, 369]}
    #kludge={'peak_id':[ 368, 369]}
    kludge={'peak_id':[367,368]}
    leaf_storage={}
    MT=mountain_top.cut_mountain_top( this_simname, target_fname = mountain_top_name, do_projections=True, verify=verify_cores_u302,kludge=kludge, leaf_storage=leaf_storage)
                                     #density_cut_fname="u302_densiy_cut_test.h5")
    overlap_tmp = check_overlap(leaf_storage)




if 'new_thresholds' not in dir():
    new_thresholds={}
if 'all_peaks' not in dir():
    all_peaks=defaultdict(dict)
skip_list = [67, 92, 102, 194, 368]

if 1:
    mountain_top.split_all(this_simname, new_thresholds, all_peaks, skip_list=skip_list)

cut_name = 'u302_density_cuts_auto2.h5'
mountain_top.dump_cutoffs(cut_name, new_thresholds)

if 1:
    kludge={'peak_id':peak_list}
    test_leaf_storage={}
    MT=mountain_top.cut_mountain_top( this_simname, target_fname = mountain_top_name, do_projections=True, verify=verify_cores_u302,kludge=kludge, leaf_storage=test_leaf_storage,
                                     density_cut_fname=cut_name)
    overlap_maybe_great = mountain_top.check_overlap(test_leaf_storage)

if 0:
    more_leaf_storage={}
    split_density=mountain_top.split_peaks(this_simname, 367,368 ,do_projections=True,leaf_storage=more_leaf_storage)
    #{17: 906.8142319419817, 19: 906.8142319419817}
    #{18: 273.2661243303226, 67: 273.2661243303226}
    # {366: 542.7954407207387, 367: 542.7954407207387}
    #{367: 75071.89940015283, 368: 75071.89940015283}

