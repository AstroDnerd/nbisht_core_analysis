
from starter2 import *
import data_locations as dl
from collections import defaultdict

import davetools
reload(davetools)

plt.close('all')


class close_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
        self.overlap=defaultdict(list)
        self.distance=defaultdict(list)
        self.last_point=[]

    def make_distance(self,core_list=None):
        dx=1./2048
        nx = 1./dx
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores
        if hasattr(core_list,'v'):
            core_list=core_list.v #needs to not have unit.
            core_list=core_list.astype('int')

        tsorted = thtr.times
        for core_id in core_list:
            
            print('distance: ', core_id)
            self.cores_used.append(core_id)
            ms = trackage.mini_scrubber(thtr,core_id)
            self.last_point.append(  ms.mean_center[:,-1])
            self.ms = ms
        self.distance_matrix = np.zeros( [len(self.cores_used)]*2) -1 #start with negative numbers for error checking
        for nc1,core_id_1 in enumerate(self.cores_used):
            for nc2,core_id_2 in enumerate(self.cores_used):
                x1 = self.last_point[nc1]
                x2 = self.last_point[nc2]
                self.distance_matrix[ nc1, nc2 ] = np.sqrt( ((x1-x2)**2).sum() )


