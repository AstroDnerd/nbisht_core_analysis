
from starter2 import *
import matplotlib.image as mpimg

import data_locations as dl
reload(dl)
reload(trackage)
plt.close('all')

form='pdf'

import means_etc
reload(means_etc)
if 'm1' not in dir():
    m1 = means_etc.means_etc(tl.looper1 )
    m2 = means_etc.means_etc(tl.looper2 )
    m3 = means_etc.means_etc(tl.looper3 )
    mmm = {'u201':m1, 'u202':m2, 'u203':m3}

   
fig, ax = plt.subplots(1,1)
colors = {'u201':'r','u202':'g','u203':'b'}
for name in mmm:
    ax.scatter( mmm[name].npart, mmm[name].dmeans, c=colors[name])
axbonk(ax,xscale='log',yscale='log',xlabel='N',ylabel='rho')
fig.savefig('plots_to_sort/density_npart.pdf')

