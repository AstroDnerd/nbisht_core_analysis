from starter2 import *
import convex_hull_tools as CHT
import three_loopers_u500 as TL
import supersets
sim_list=['u501','u502','u503']

if 'ht' not in dir() :
    ht = {}
    for this_simname in sim_list:
        ht[this_simname] = CHT.hull_tool(TL.loops[this_simname])
        ht[this_simname].make_hulls()
        ht[this_simname].make_overlaps()
if 'st' not in dir():
    st={}
    for this_simname in sim_list:
        st[this_simname] = supersets.superset( TL.loops[this_simname], ht[this_simname])
        st[this_simname].find()

for sim in sim_list:
    this_looper = TL.loops[sim]
    core_list=np.unique(this_looper.tr.core_ids)
    core_id_list=[]
    value_list=[]
    for core_id in sorted(core_list):
        core_id_list.append(core_id)
        value_list.append( st[sim].set_by_core[core_id])
    fptr = h5py.File('browser_data/neighborhood_by_core_%s.h5'%sim,'w')
    fptr['core_ids']=nar(core_id_list)
    fptr['neighborhood']=nar(value_list)
    fptr.close()
