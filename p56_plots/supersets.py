
from starter2 import *

import three_loopers_mountain_top as TLM
import convex_hull_tools as CHT

if 'clobber' not in dir():
    clobber=False

def recursive_add_cores(core_id, superset=None, thresh=0.0, cores_used=None, overlap_matrix=None):
    if superset is None:
        superset = set()
    superset.add(core_id)

    more_index = np.where( cores_used == core_id )[0][0]
    #more_overlaps = (overlap_matrix[more_index,:].flatten() >thresh) + (overlap_matrix[:,more_index].flatten()>thresh)
    both_directions = np.column_stack([overlap_matrix[more_index,:].flatten(), overlap_matrix[:,more_index].flatten()])
    this_overlap = both_directions.min(axis=1)
    more_overlaps = this_overlap > thresh
    more_cores = cores_used[ more_overlaps]
    for core in more_cores:
        if core not in superset:
            recursive_add_cores(core,  superset=superset, cores_used=cores_used, overlap_matrix=overlap_matrix)
    return superset

class superset():
    def __init__(self,this_looper,hull ):
        self.this_looper=this_looper
        self.name = this_looper.sim_name
        self.supersets = []
        self.set_by_core = {}
        self.hull = hull
    def find(self, thresh=0.0):

        #values from cores_all get popped out and used in supersets.
        cores_all = set( self.hull.cores_used)
        cores_used = nar(self.hull.cores_used)

        overlap_matrix = np.zeros( [len(cores_used)]*2)
        for ncore, core_id in enumerate(cores_used):
            overlap_matrix[ncore,:] = nar(self.hull.overlaps[core_id])

        while len(cores_all):
            this_core = cores_all.pop()
            print("===",this_core,"===")
            this_superset = recursive_add_cores(this_core, overlap_matrix=overlap_matrix, cores_used=cores_used, thresh=thresh)
            print(this_superset)
            self.supersets.append(this_superset)
            cores_all.difference_update(this_superset)
        for nS, superset in enumerate(self.supersets):
            for core_id in superset:
                self.set_by_core[core_id] = nS
if 0:
    for nset,this_superset in enumerate(supersets):
        CHT.plot_2d(HT,core_list=this_superset,frames=[0],accumulate=True,label_cores=[0], prefix = "S%02d"%nset)

