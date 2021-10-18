from starter2 import *

import mountain_top
reload(mountain_top)
import three_loopers_mountain_top as TLM
frame=0
ds = TLM.loops['u301'].load(frame)
trees = mountain_top.target_tree(ds,"datasets_small/u301_mountain_tops_take_9.h5")
trees.make_leaves_check()
