from starter2 import *

import three_loopers_u500 as TL

sim_list = ['u501','u502','u503']
sim_list=['u501']

centers={}
zooms={}
centers['u503'] = [0.6,0.78]
zooms['u503'] = 4
centers['u502']=[0.1,0.2]
zooms['u502']=3
centers['u501']=[0.425,0.475]
zooms['u501']=5
for sim in sim_list:
    if sim not in centers:
        continue
    loop = TL.loops[sim]
    ds = loop.load( loop.target_frame )
    pw = yt.ProjectionPlot( ds, 0, 'density', origin='domain')
    pw.set_center( centers[sim])
    pw.zoom( zooms[sim])
    pw.save('plots_to_sort/thumbnail_%s'%sim)
