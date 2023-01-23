
from starter2 import *

import three_loopers_u500 as TL


for sim in TL.loops:
    this_looper=TL.loops[sim]
    print( this_looper.core_ids.size, this_looper.core_by_mode['One'].size, this_looper.core_by_mode['Binary'].size, this_looper.core_by_mode['Cluster'].size)
    print( this_looper.core_ids.size, this_looper.core_by_mode['One'].size+ this_looper.core_by_mode['Binary'].size+ this_looper.core_by_mode['Cluster'].size)

