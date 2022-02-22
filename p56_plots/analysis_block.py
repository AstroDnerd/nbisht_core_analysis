

from starter2 import *

import convex_hull_tools as CHT
import matplotlib.colors as colors

reload(CHT)
import hair_dryer
reload(hair_dryer)
import stay_close
import three_loopers_tenfour as TL4
sim_list=['u401','u402','u403']

import three_loopers_tenfour as TL4
import convex_hull_tools as CHT
reload(CHT)
sim_list=['u401','u402','u403']
if 'ht' not in dir() :
    ht = {}
    for this_simname in sim_list:
        ht[this_simname] = CHT.hull_tool(TL4.loops[this_simname])
        ht[this_simname].make_hulls()
        ht[this_simname].make_overlaps()
if 'ct' not in dir():
    ct = {}
    for this_simname in sim_list:
        ct[this_simname] = stay_close.close_tool( TL4.loops[this_simname])
        ct[this_simname].make_distance()

import supersets
reload(supersets)
if 'st' not in dir():
    st={}
    for this_simname in sim_list:
        st[this_simname] = supersets.superset( TL4.loops[this_simname], ht[this_simname])
        st[this_simname].find()
