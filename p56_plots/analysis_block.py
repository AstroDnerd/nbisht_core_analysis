"""
Analysis block 1
"""

from starter2 import *

import convex_hull_tools as CHT
import matplotlib.colors as colors

import r_inflection
reload(r_inflection)
import supersets
reload(supersets)
reload(CHT)
import hair_dryer
reload(hair_dryer)
import close_tool
import three_loopers_six as TL
sim_list=['u601','u602','u603']

if 'ht' not in dir():
    ht = {}
    for this_simname in sim_list:
        ht[this_simname] = CHT.hull_tool(TL.loops[this_simname])
        ht[this_simname].make_hulls()
        ht[this_simname].make_overlaps()
if 'ct' not in dir():
    ct = {}
    for this_simname in sim_list:
        ct[this_simname] = close_tool.close_tool( TL.loops[this_simname])
        ct[this_simname].make_distance()

if 'st' not in dir():
    st={}
    for this_simname in sim_list:
        st[this_simname] = supersets.superset( TL.loops[this_simname], ht[this_simname])
        st[this_simname].find()

if 'inflection' not in dir():
    inflection = {}
for sim in sim_list:
    if sim in inflection:
        continue
    inflection[sim]=r_inflection.R_INFLECTION( TL.loops[sim])
    inflection[sim].run()

