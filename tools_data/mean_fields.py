

from starter2 import *
import three_loopers_u500 as TL

means={}
for sim in TL.loops:
    this_looper=TL.loops[sim]
    ds = this_looper.load(0)
    ad = ds.all_data()
    bx = ad['magnetic_field_x'].mean()
    by = ad['magnetic_field_y'].mean()
    bz = ad['magnetic_field_z'].mean()
    means[sim]=[bx,by,bz]
    print(sim,means[sim])
    

