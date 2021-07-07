
from starter2 import *
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
reload(dl)
plt.close('all')
color={'u05':'r','u10':'g','u11':'b'}
color.update({'u201':'r','u202':'g','u203':'b'})

import convex_hull_tools as CHT
reload(CHT)

import three_loopers_1tff as tl

if 'clobber' not in dir():
    clobber=False

    
if 'ht1' not in dir() or clobber: 
    ht1 = CHT.hull_tool(tl.looper1)
if 'ht2' not in dir() or clobber:
    ht2 = CHT.hull_tool(tl.looper2)
    #ht2.plot_2d(frames=[0])
if 'ht3' not in dir() or clobber:
    ht3 = CHT.hull_tool(tl.looper3)

for frame in [0]:
    if 0:
        #
        # Compute overlaps
        #
        for htool in [ht3]: #[ht1, ht2, ht3]:
            htool.overlaps=defaultdict(list)
            htool.make_hulls(frames=[frame])
            for core_1 in htool.cores_used:
                print("overlap li,", core_1)
                for core_2 in htool.cores_used:
                    result = htool.check_hull_overlap(core_1,core_2)
                    if core_1 == core_2:
                        result = -result
                    htool.overlaps[core_1].append(result)
        fractions,cores_84=CHT.get_overlapping_cores(ht3,84)


if 1:
    for frame in [0] :#range(50,110,10):
        if 0:
            #
            # Compute hulls
            #
            for htool in [ht3]: #[ht1, ht2, ht3]:
                htool.overlaps=defaultdict(list)
                htool.make_hulls(frames=[frame])
#84, 173, 210
        if 1:
            ht3b = CHT.hull_tool(tl.looper3)
            CHT.plot_2d(ht3b,core_list = cores_84, accumulate=True,frames=[frame], all_plots=True)

#
# core time series
#
if 0:
    for frame in range(0,110,10):
        ht3b = CHT.hull_tool(tl.looper3)
        ht3b.make_hulls( core_list = [128])
        CHT.plot_2d(ht3b,core_list = [128], accumulate=True,frames=[frame])

