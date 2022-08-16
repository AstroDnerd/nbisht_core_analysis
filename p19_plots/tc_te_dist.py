
from starter2 import *
from collections import defaultdict
import scipy
import colors

from scipy.ndimage import gaussian_filter
import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL
import movie_frames 
from collections import defaultdict

import tsing
reload(tsing)
def cuml(ax,quan,color,label):
    the_x = nar(sorted(quan))
    the_y = np.arange( the_x.size)/the_x.size
    ax.plot( the_x, the_y, color=color,label=label)
sims=['u501', 'u502','u503']
#sims=[ 'u502','u503']
#sims=['u501']
for sim in sims:
    core_list=None
if 1:
    objects={}
    fig,axes=plt.subplots(3,3)
    fig.subplots_adjust(wspace=0,hspace=0)
    fsung_ext=extents()
    fsung_d={}
    for ns,sim in enumerate(sims):
        obj=tsing.te_tc(TL.loops[sim])
        objects[sim]=obj
        obj.run()
        for nq,Q in enumerate(['One','Binary','Cluster']):
            color = 'rgb'[nq]
            cuml( axes[ns][0], obj.tsing[Q], color=color,label=Q)
            cuml( axes[ns][1], nar(obj.tend[Q])-nar(obj.tsing[Q]), color=color, label=Q)
            #axes[ns][0].hist( obj.tsing[Q], color=color, histtype='step', label=Q)
            #axes[ns][1].hist( , color=color, histtype='step', label=Q)
        axes[ns][0].set(xlabel=r'$t_{\rm{sing}}$',xlim=[0,0.9], ylabel='N sim %s'%(ns+1), ylim=[0,1], yticks=[0.25,0.5,0.75], xticks=[0.25,0.5,0.75])
        axes[ns][1].set(xlabel=r'$\Delta t_{\rm{sing}}$',xlim=[0.01,0.2], ylim=[0,1],yticks=[], xticks=[0.05,0.1,0.15])
        fsung = nar(obj.fsung)
        fsung_ext( fsung[fsung>0])
        fsung_d[sim]=fsung



        #plt.
        #plt.hist( obj.tsing['Binary'], color='g', histtype='step', label='Binary')
        #plt.hist( obj.tsing['Cluster'], color='b', histtype='step', label='Cluster')

    for ns,sim in enumerate(sims):
        obj = objects[sim]
        for nq,Q in enumerate(['One','Binary','Cluster']):
            color = 'rgb'[nq]
            mask = nar(obj.mode) == (nq+1)
            bins = np.geomspace( fsung_ext.minmax[0]*(1+1*nq), 1, 16)
            bins = np.linspace( -0.01*nq, 1, 16)
            fsung = nar(obj.fsung)[mask]
            cuml(axes[ns][2],fsung, color=color, label=None)

            #axes[ns][2].hist( fsung[fsung>0], bins=bins,histtype='step', color=color)
        axes[ns][2].set(xscale='linear', yticks=[], xlabel=r'$f_{sing}$', xticks=[0.1,0.3,0.5, 0.7, 0.9])
        if ns<2:
            axes[ns][0].set(xticks=[])
            axes[ns][1].set(xticks=[])
            axes[ns][2].set(xticks=[])
        if ns==0:
            axes[ns][0].legend(loc=2)
    fig.savefig('plots_to_sort/tsing.pdf')

