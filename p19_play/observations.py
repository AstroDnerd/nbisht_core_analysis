
from starter2 import *

import pgcc_tools as pgcc
reload(pgcc)

types={'name':str, 'ra':str, 'dec':str,'temperature':float,'temperature_error':float, 'mass':float,
       'mass_error':float,'size':float,'size_error':float}
stuff=pgcc.parse_results("p19_plots/PGCC_raw.txt",types)
#print(stuff.data.keys())
mask=stuff.return_multi_mask(['mass','size'])

fig,axes=plt.subplots(2,2)
ax0=axes[0][0];ax1=axes[0][1]
ax2=axes[1][0]

masses = nar(stuff.data['mass'])[mask]
sizes  = nar(stuff.data['size'])[mask]
mbins = np.geomspace(masses.min(),masses.max(),10)
rbins = np.geomspace(sizes.min(),sizes.max(),10)
mask2 = (sizes<20)
ax0.hist( masses, bins=mbins)
ax1.hist( sizes, bins=rbins)
ax0.set(xscale='log')
ax1.set(xscale='log')
ax2.scatter( sizes, masses)
ax2.set(xscale='log',yscale='log')

fig.savefig('%s/hist.pdf'%plot_dir)

