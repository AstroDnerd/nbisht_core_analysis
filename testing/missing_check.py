from starter2 import *


import three_loopers_mountain_top as TLM

for this_simname in TLM.loops:
    density = TLM.loops[this_simname].tr.track_dict['density']
    if density.min() <= 0:
        print("BAD ERROR sim %s has missing density"%this_simname)

