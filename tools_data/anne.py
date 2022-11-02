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
def make_hulls():
    for this_simname in sim_list:
        if sim in ht:
            continue
        print("Analysis: Hulls")
        ht[this_simname] = CHT.hull_tool(TL.loops[this_simname])
        ht[this_simname].make_hulls()
        ht[this_simname].make_overlaps()

if 'ct' not in dir():
            ct = {}
def make_close_tool():
    for this_simname in sim_list:
        if sim in ct:
            continue
        print("Analysis, close")
        ct[this_simname] = close_tool.close_tool( TL.loops[this_simname])
        ct[this_simname].make_distance()

if 'st' not in dir():
    st={}
def make_subsets():
    for this_simname in sim_list:
        if sim in st:
            continue
        print("Analysis: supersets")
        st[this_simname] = supersets.superset( TL.loops[this_simname], ht[this_simname])
        st[this_simname].find()

if 'inflection' not in dir():
    inflection = {}
def make_inflection():
    for sim in sim_list:
        if sim in inflection:
            continue
        print("Analysis: inflection")
        inflection[sim]=r_inflection.R_INFLECTION( TL.loops[sim])
        inflection[sim].run()

def dump_inflection():
    for sim in inflection:
        fname = "browser_data/r_inflection_%s.h5"%sim
        if len(glob.glob(fname)) > 0:
            print("Will not clobber %s"%fname)
            continue
        core_ids = sorted(inflection[sim].rinflection.keys())
        rrr = [inflection[sim].rinflection[core_id] for core_id in core_ids]
        fptr = h5py.File(fname,'w')
        fptr['core_ids'] = core_ids
        fptr['r_inflection'] = rrr
        fptr.close()

class r_infl():
    #a fake inflection class.
    def __init__(self, core_ids, r_inflection):
        self.rinflection = dict( zip( core_ids, r_inflection))
def read_inflection():
    for sim in sim_list:
        fname = "browser_data/r_inflection_%s.h5"%sim
        if len(glob.glob(fname)) == 0:
            print("Can not clobber %s"%fname)
            continue
        fptr = h5py.File(fname,'r')
        core_ids = fptr['core_ids'][()]
        r_inflection = fptr['r_inflection'][()]
        fptr.close()
        inflection[sim] = r_infl(core_ids, r_inflection)

