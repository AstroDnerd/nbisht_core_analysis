
'''
three way plots to see how alpha 1 vs alpha 2 relate...
'''
from starter2 import *
import means_etc

figa,axa,axtop,axright = means_etc.three_way_bean()
axa.scatter( NumberOfOverlap, Fraction,c='k')  #I don't think we need this 'b'
axa.scatter( NumberOfOverlap,Max ,c='r')

axtop.hist( NumberOfOverlap, bins=bins_n, histtype='step',color='k')
axright.hist( Fraction, bins=bins_f, histtype='step',color='k',orientation='horizontal')
axright.hist( Max, bins=bins_f, histtype='step',color='r',orientation='horizontal')

#QUESTION: what was the purpose of plotting MAX here?


if 1:
    axbonk(axa,xlabel=r'$N_{\rm{overlap}}$', ylabel='Overlap Fraction')
    axa.set_xlim([-0.1,nmax+0.1])
    axbonk(axtop,xlabel='',ylabel=r'$N$')
    axbonk(axright,xlabel=r'$N$',ylabel='')
    axright.set_ylim( axa.get_ylim())
    axtop.set_xlim( axa.get_xlim())
    axtop.set_xticks([])
    axright.set_yticks([])

