
# three_way plotting
from starter2 import *
import means_etc
reload(means_etc)

plt.clf()

figa, axa, axtop, axright = means_etc.three_way_bean()

axa.scatter( NumberOfOverlap, Fraction,c='k')
axtop.hist( NumberOfOverlap, bins=bins_n, histtype='step',color='k')
axright.hist( Fraction, bins=bins_f, histtype='step',color='k',orientation='horizontal')

axa.set_xlim([-0.1,nmax+0.1])
axright.set_ylim(axa.get_ylim())
axtop.set_xlim(axa.get_xlim())
axtop.set_xticks([])
axright.set_yticks([])



