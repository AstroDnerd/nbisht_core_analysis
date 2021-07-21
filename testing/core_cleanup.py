from starter2 import *
import mountain_top
reload(mountain_top)

#run core_multi first.
#needs overlap_maybe_great and overlap_all to be populated.

all_to_check={}
for peak_id in overlap_maybe_great:
    if peak_id != 43:
        continue
    all_to_check[peak_id] = overlap_maybe_great[peak_id]
    all_to_check[peak_id] +=overlap_all[peak_id]

if 'new_thresholds_cleanup' not in dir() or True:
    new_thresholds_cleanup={}
if 'all_peaks_cleanup' not in dir() or True:
    all_peaks_cleanup=defaultdict(dict)


if 1:
    mountain_top.split_all(this_simname, all_to_check, new_thresholds_cleanup, all_peaks_cleanup, do_projections=True)

temp_threshold_name = "u302_cutoff_cleanup.h5"
os.remove(temp_threshold_name)
mountain_top.dump_cutoffs(temp_threshold_name,new_thresholds_cleanup)

if 1:
    print("check the density in the file")
    pdb.set_trace()
    peak_list = list(all_to_check.keys())
    kludge={'peak_id':peak_list}
    clean_leaf_storage={}
    temp_target = 'u302_mountain_tops_cleanupb.h5'
    MT=mountain_top.cut_mountain_top( this_simname, target_fname = temp_target, do_projections=True, verify=verify_cores_u302,kludge=kludge, leaf_storage=clean_leaf_storage,
                                     density_cut_fname=temp_threshold_name)
    overlap_cleanup = mountain_top.check_overlap(clean_leaf_storage)

