
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



import three_loopers_tenfour as TL4
import convex_hull_tools as CHT
reload(CHT)
sim_list=['u401']#,'u402','u403']
if 'ht' not in dir() :
    ht = {}
    for this_simname in sim_list:
        ht[this_simname] = CHT.hull_tool(TL4.loops[this_simname])
        ht[this_simname].make_hulls()
        ht[this_simname].make_overlaps()
if 'ct' not in dir():
    ct = {}
    for this_simname in sim_list:
        ct[this_simname] = close_tool( TL4.loops[this_simname])
        ct[this_simname].make_distance()

import supersets
reload(supersets)
if 'st' not in dir():
    st={}
    for this_simname in sim_list:
        st[this_simname] = supersets.superset( TL4.loops[this_simname], ht[this_simname])
        st[this_simname].find()

if 1:
    #every pair
    import means_etc
    for this_simname in sim_list:
        htool = ht[this_simname]
        ctool = ct[this_simname]
        stool = st[this_simname]
        overlap_matrix = np.zeros( [len(htool.cores_used)]*2) -1
        overlap_matrix_anti = np.zeros( [len(htool.cores_used)]*2) -1
        distance_matrix = np.zeros( [len(htool.cores_used)]*2) -1
        color_matrix = np.zeros_like(distance_matrix).astype('int')
        for nc1,core_id_1 in enumerate(htool.cores_used):
            for nc2,core_id_2 in enumerate(htool.cores_used):
                #if nc2 <= nc1:
                #    continue
                overlap_matrix_anti[nc1,nc2] = htool.overlaps[core_id_1][nc2]
                overlap_matrix[nc1,nc2] = min( [htool.overlaps[core_id_1][nc2], htool.overlaps[core_id_2][nc1] ])
                color = 1
                if stool.set_by_core[core_id_1] == stool.set_by_core[core_id_2]:
                    color = 2
                color_matrix[nc1,nc2]=color
        #fig,ax=plt.subplots(1,1)
        fig, ax, axtop,axright = means_etc.three_way_bean()
        Md = ctool.distance_matrix.flatten()
        Mo = overlap_matrix.flatten()
        ok = (Md >= 0)*(Mo >= 0)
        color = (color_matrix).flatten()[ok]
        c = ['grk'[v] for v in color]
        ax.scatter( Md[ok][ color==1], Mo[ok][ color==1]-0.01,color='r')
        ax.scatter( Md[ok][ color==2], Mo[ok][ color==2],color='k')
        bins = np.geomspace( Md[ok].min(), Md[ok].max(),64)
        axtop.hist( Md[ok][ color == 1], histtype='step',color='r',bins=bins)
        axtop.hist( Md[ok][ color == 2], histtype='step',color='k',bins=bins)
        #axright.hist( Mo[ok][color == 1], histtype='step',color='r',orientation='horizontal',bins=64)
        axright.hist( Mo[ok][color == 2], histtype='step',color='k',orientation='horizontal',bins=64)
        axbonk(ax, xlabel='Distance', ylabel='overlap',xscale='log',yscale='linear')
        axbonk(axtop, xlabel=None, ylabel='N',xscale='log',yscale='log')
        axbonk(axright, xlabel='N', ylabel=None,xscale='log',yscale='linear')
        axtop.set_xlim( ax.get_xlim())
        axright.set_ylim( ax.get_ylim())

        fig.savefig('plots_to_sort/%s_stay_close.png'%(htool.name))
        plt.close('all')


if 0:
    #closest pairs
    import means_etc
    for this_simname in sim_list:
        htool = ht[this_simname]
        ctool = ct[this_simname]
        stool = st[this_simname]
        overlap_matrix = np.zeros( [len(htool.cores_used)]*2) -1
        distance_matrix = np.zeros( [len(htool.cores_used)]*2) -1
        color_matrix = np.zeros_like(distance_matrix).astype('int')
        for nc1,core_id_1 in enumerate(htool.cores_used):
            for nc2,core_id_2 in enumerate(htool.cores_used):
                overlap_matrix[nc1,nc2] = max( [htool.overlaps[core_id_1][nc2], htool.overlaps[core_id_2][nc1] ])
                color = 1
                if stool.set_by_core[core_id_1] == stool.set_by_core[core_id_2]:
                    color = 2
                color_matrix[nc1,nc2]=color
        #fig,ax=plt.subplots(1,1)
        fig, ax, axtop,axright = means_etc.three_way_bean()
        Md1 = ctool.distance_matrix +0
        Mo1 = overlap_matrix

        Md1[ Md1==0] = 5 #the closest core to any one is itself.  So make zero 5 so we can take the min.
        closest = np.argmin(Md1,axis=0)

        Md = Md1[closest, tuple(range(Md1.shape[1]))]
        Mo = Mo1[closest, tuple(range(Md1.shape[1]))]

if 0:

        #Md = ctool.distance_matrix.flatten()
        #Mo = overlap_matrix.flatten()
        #ok = (Md >= 0)*(Mo >= 0)
        ok = slice(None)
        color = (color_matrix[closest, tuple(range(Md1.shape[1]))]).flatten()[ok]
        c = ['grk'[v] for v in color]
        ax.scatter( Md[ok][ color==1], Mo[ok][ color==1]-0.01,color='r')
        ax.scatter( Md[ok][ color==2], Mo[ok][ color==2],color='k')
        bins = np.geomspace( Md[ok].min(), Md[ok].max(),64)
        axtop.hist( Md[ok][ color == 1], histtype='step',color='r',bins=bins)
        axtop.hist( Md[ok][ color == 2], histtype='step',color='k',bins=bins)
        #axright.hist( Mo[ok][color == 1], histtype='step',color='r',orientation='horizontal',bins=64)
        axright.hist( Mo[ok][color == 2], histtype='step',color='k',orientation='horizontal',bins=64)
        axbonk(ax, xlabel='Distance', ylabel='overlap',xscale='log',yscale='linear')
        axbonk(axtop, xlabel=None, ylabel='N',xscale='log',yscale='log')
        axbonk(axright, xlabel='N', ylabel=None,xscale='log',yscale='linear')
        axtop.set_xlim( ax.get_xlim())
        axright.set_ylim( ax.get_ylim())

        fig.savefig('plots_to_sort/%s_stay_closest.png'%(htool.name))
        plt.close('all')






