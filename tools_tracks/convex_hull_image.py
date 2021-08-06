
from starter2 import *
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
reload(dl)
plt.close('all')
color={'u05':'r','u10':'g','u11':'b'}
color.update({'u201':'r','u202':'g','u203':'b'})
color.update({'u301':'r','u302':'g','u303':'b'})

import convex_hull_tools as CHT
reload(CHT)

#import three_loopers_1tff as tl
import three_loopers_mountain_top as TLM

if 'clobber' not in dir():
    clobber=False

if 'ht1' not in dir() or clobber: 
    ht1 = CHT.hull_tool(TLM.loops['u301'])
if 'ht2' not in dir() or clobber:
    ht2 = CHT.hull_tool(TLM.loops['u302'])
    #ht2.plot_2d(frames=[0])
if 'ht3' not in dir() or clobber:
    ht3 = CHT.hull_tool(TLM.loops['u303'])

#test versions
if 0:
    CHT.plot_2d(ht1,frames=[0],core_list=[85,86, 306, 307, 308], accumulate=True)
if 0:
    CHT.plot_2d(ht1,core_list=[24],frames=[0], accumulate=True, label_cores=[-1])
if 0:
    CHT.plot_2d(ht2,frames=[0], accumulate=True, label_cores=[-1])
if 0:
    CHT.plot_2d(ht3,frames=[0], accumulate=True, label_cores=[])

#paper versions
if 0:
    CHT.plot_2d(ht1,frames=[0], accumulate=True, label_cores=[323])
if 1:
    CHT.plot_2d(ht2,frames=[0], accumulate=True, label_cores=[])
if 1:
    CHT.plot_2d(ht3,frames=[0], accumulate=True, label_cores=[])
