
from starter2 import *
from collections import defaultdict
import scipy
import colors

from scipy.ndimage import gaussian_filter
import hair_dryer
reload(hair_dryer)

import track_loader as TL
import movie_frames 
from collections import defaultdict

import tsing
reload(tsing)
def cuml(ax,quan,color,label):
    the_x = nar(sorted(quan))
    the_y = np.arange( the_x.size)/the_x.size
    ax.plot( the_x, the_y, color=color,label=label)
sims=['u501', 'u502','u503']
track_loader.load_tracks(sims)
#sims=[ 'u502','u503']
#sims=['u501']
#sims=['u502']
for sim in sims:
    core_list=None
if 1:
    objects={}
    fig,axes=plt.subplots(3,2)
    fig.subplots_adjust(wspace=0,hspace=0)
    fsung_ext=extents()
    fsung_d={}
    for ns,sim in enumerate(sims):
        obj=tsing.te_tc(TL.loops[sim])
        objects[sim]=obj
        obj.run()
        for nq,Q in enumerate(['Alone','Binary','Cluster']):
            color = 'rgb'[nq]
            cuml( axes[ns][0], obj.tsing[Q], color=color,label=Q)
            cuml( axes[ns][1], nar(obj.tend[Q])-nar(obj.tsing[Q]), color=color, label=Q)
            #axes[ns][0].hist( obj.tsing[Q], color=color, histtype='step', label=Q)
            #axes[ns][1].hist( , color=color, histtype='step', label=Q)
        axes[ns][0].set(xlabel=r'$t_{\rm{sing}}/t_{\rm{ff}}$',xlim=[0,0.9], ylabel=r'$N/N_{total}$', ylim=[0,1], yticks=[0.25,0.5,0.75], xticks=[0.25,0.5,0.75])
        axes[ns][1].set(xlabel=r'$\Delta t_{\rm{sing}}/t_{\rm{ff}}$',xlim=[0.01,0.2], ylim=[0,1],yticks=[], xticks=[0.05,0.1,0.15])
        fsung = nar(obj.fsung)
        fsung_ext( fsung[fsung>0])
        fsung_d[sim]=fsung

    for ns,sim in enumerate(sims):
        obj = objects[sim]
        axes[ns][0].text(0.1,0.75,r'$sim %d$'%(ns+1))
        if ns<2:
            axes[ns][0].set(xticks=[])
            axes[ns][1].set(xticks=[])
            #axes[ns][2].set(xticks=[])
        if ns==2:
            axes[ns][1].legend(loc=4)
    fig.savefig('plots_to_sort/tsing.pdf')

