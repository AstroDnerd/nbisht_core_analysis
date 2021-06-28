from starter2 import *

import three_loopers_1tff as tl
if 0:
    for loop in [tl.looper1, tl.looper2, tl.looper3]:
        for frame in loop.snaps:
            for core_id in loop.snaps[frame]:
                loop.snaps[frame][core_id].ind = None
                loop.snaps[frame][core_id].vel = None
        output_name = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/every_ten_u20x_compact/%s.h5"%loop.out_prefix
        loop.save(output_name)

import three_loopers_take3
reload(three_loopers_take3)
