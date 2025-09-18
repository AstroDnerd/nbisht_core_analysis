
from starter2 import *
from collections import defaultdict
import heat_map
from scipy.optimize import curve_fit
import colors
reload(heat_map)

class vheat():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
        self.giant_array=None

    def run(self,core_list=None):
        thtr=self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        self.times=tsorted

        self.giant_array=np.zeros_like(thtr.track_dict['velocity_x'])
        self.vx_array=np.zeros_like(thtr.track_dict['velocity_x'])
        npart=0
        for core_id in core_list:
            self.cores_used.append(core_id)
            ms = trackage.mini_scrubber(thtr,core_id)
            nthis = ms.nparticles
            self.giant_array[npart:npart+nthis,:] = ms.rel_vmag
            self.vx_array[npart:npart+nthis,:] = ms.rel_vx
            npart = npart + nthis


import three_loopers_u500 as TL

simlist = ['u501','u502','u503']
simlist = ['u501']

if 'vtools' not in dir():
    vtools={}
    for sim in simlist:
        vtools[sim]=vheat( TL.loops[sim])
        vtools[sim].run()

if 1:
    fig,axes=plt.subplots(1,3, figsize=(12,8), sharey=True)
    axlist=axes.flatten()
    for ns,sim in enumerate(vtools):
        ax=axlist[ns]
        this_looper=vtools[sim].this_looper
        Q = vtools[sim].vx_array
        mx = Q.max()
        mn = Q[Q>1e-8].min()
        minmin = min([mn, 1/mx])
        maxmax = max([mx, 1/mn])

        bins = np.geomspace( minmin,maxmax,64)
        heat_map.heat_map(Q[:,:-1], this_looper.tr.times[:-1], ax=ax,bins=bins)
        axbonk(ax,yscale='log',ylim=[minmin,maxmax])
    fig.savefig('plots_to_sort/velocity_heat.png')

