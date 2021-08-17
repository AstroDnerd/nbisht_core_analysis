from starter2 import *

#import three_loopers_mountain_top as TLM
#
#for this_simname in TLM.loops:
#    print(this_simname, (TLM.loops[this_simname].tr.track_dict['density'] == 0).any())

for core_id in new_looper.snaps[50]:
    snap5 = new_looper.snaps[50][core_id]
    snap4 = new_looper.snaps[40][core_id]
    index_error=np.abs(snap5.ind - snap4.ind)
    if index_error.any():
        print("OH NO")
