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
    skip += [67] #diffuse, hard to separate from 66
    if peak_id in skip:
        out=False
    return out

#
# Make all the leaves and don't touch it.
#
#if 'leaf_storage_all' not in dir():#'leaf_storage' not in dir() and False:
#    kludge={}
#    leaf_storage_all={}
#    MT=mountain_top.cut_mountain_top( this_simname, target_fname = mountain_top_name, do_projections=do_mountain_projections, verify=verify_cores_u302,kludge=kludge, leaf_storage=leaf_storage_all)
#                                     #density_cut_fname="u302_contour_mins.h5")
#    overlap_all=check_overlap(leaf_storage_all)

    

if 1:
    outname = "u302_tracks_core_258.h5"
#kludge={'peak_id':[0,1,364, 113,u302_was258]})
#kludge={'peak_id':[258]}
#kludge={'peak_id':[core_id]}
    kludge={}
#kludge={'peak_id':[89,11,14,350, 349, 368, 369]}
#kludge={'peak_id':[ 368, 369]}
    k0=186
    k1=188   
    kludge={'peak_id':[k0,k1]}
    leaf_storage_here={}
    MT=mountain_top.cut_mountain_top( this_simname, target_fname = mountain_top_name, do_projections=True, verify=verify_cores_u302,
                                     kludge=kludge, leaf_storage=leaf_storage_here, radius_dict=coreset_data.radius_302)
                                     #density_cut_fname="u302_densiy_cut_working.h5")
    overlap_tmp = mountain_top.check_overlap(leaf_storage_here)

reload(coreset_data)

if 1:
    more_leaf_storage_here={}
    split_density=mountain_top.split_peaks(this_simname, k0,k1,do_projections=True,
                                           leaf_storage=more_leaf_storage_here, radius_dict=coreset_data.radius_302,recheck=True)
    #{17: 906.8142319419817, 19: 906.8142319419817}
    #{18: 273.2661243303226, 67: 273.2661243303226}
    # {366: 542.7954407207387, 367: 542.7954407207387}
    #{367: 75071.89940015283, 368: 75071.89940015283}


if 0:
    outname = "TEMP6.h5"
    kludge={}
    kludge={'peak_id':[k0,k1]}
    leaf_storage_check={}
    MT=mountain_top.cut_mountain_top( this_simname, target_fname = outname, do_projections=True, verify=verify_cores_u302,
                                     kludge=kludge, leaf_storage=leaf_storage_check, cut_override=split_density)
                                     #density_cut_fname="u302_densiy_cut_working.h5")
    overlap_tmp = mountain_top.check_overlap(leaf_storage_check)




if 0:
    #it is found that 11 overlaps with 14 and 89.
    more_leaf_storage={}
    split_density=mountain_top.split_peaks(this_simname, 11,14 ,do_projections=True,leaf_storage=more_leaf_storage)
    # {11: 33027.859145733135, 14: 33027.859145733135}
density_cut_fname = "u302_density_cut_c.h5"
if 0:    
    dump_cutoffs(density_cut_fname,split_density)
    #check_overlap(more_leaf_storage)

if 0:
    #it is found that 11 overlaps with 14 and 89.
    #{18: 604.5428212946545, 17: 604.5428212946545}
    more_leaf_storage={}
    split_density=mountain_top.split_peaks(this_simname, 18,17   ,do_projections=True,leaf_storage=more_leaf_storage, radius=2e-2)
    #dump_cutoffs(density_cut_fname,split_density)

if 0:
    #3 and 5 are not overlapping, not sure why overlap_all has them.
    more_leaf_storage={}
    split_density=mountain_top.split_peaks(this_simname, 3,5 ,do_projections=True,leaf_storage=more_leaf_storage)

if 0:
    sets_so_far={18: 604.5428212946545, 17: 604.5428212946545, 11: 33027.859145733135, 14: 33027.859145733135}
    dump_cutoffs(density_cut_fname,sets_so_far)


if 0:# 'leaf_storage_sep' not in dir():#'leaf_storage' not in dir() and False:
    #kludge={'peak_id':[0,1,364, 113,u302_was258]})
    #kludge={'peak_id':[258]}
    #kludge={'peak_id':[core_id]}
    kludge={'peak_id':[11,14,17,18]}
    #kludge={'peak_id':[89,11,14,350, 349, 368, 369]}
    #kludge={'peak_id':[ 368, 369]}
    target_fname = "u302_mountain_tops_take4.h5"
    leaf_storage_sep={}
    MT=mountain_top.cut_mountain_top( this_simname, target_fname = target_fname, do_projections=True, verify=verify_cores_u302,kludge=kludge, leaf_storage=leaf_storage_sep,
                                     density_cut_fname=density_cut_fname)
    check_overlap(leaf_storage_sep)

if 0:
    leaf_store_3={}
    more_things=mountain_top.find_biggest_contour(this_simname, 11,14 , more_leaf_storage['temp_clump'],do_projections=True,leaf_storage=leaf_store_3)

if 0:
    dump_cutoffs("u302_contour_mins.h5",split_density)

