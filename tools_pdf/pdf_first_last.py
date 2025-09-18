
from starter2 import *
import data_locations as dl


def make_prof(ds,fields,weight_field=None,accumulation=False,fractional=True,n_bins=64,extrema=None):
    reg = ds.all_data()
    prof = yt.create_profile(reg,fields[0],fields[1] ,weight_field=weight_field,accumulation=accumulation,
                            fractional=fractional, n_bins=n_bins, extrema=extrema)
    the_x = 0.5*(prof.x_bins[1:]+prof.x_bins[0:-1])
    the_y = prof[fields[1]]
    #if units[0] is not None:
    #    the_x = the_x.in_units(units[0])
    #if units[1] is not None and fractional is not True:
    #    the_y = the_y.in_units(units[1])
    output={}
    output['prof']=prof
    output['the_x']=the_x
    output['the_y']=the_y
    return output


if 'prof0' not in dir():
    ds0 = yt.load("%s/DD%04d/data%04d"%(dl.enzo_directory,0,0))
    ds1 = yt.load("%s/DD%04d/data%04d"%(dl.enzo_directory,125,125))
    prof0=make_prof(ds0,['density','cell_volume'])
    prof1=make_prof(ds1,['density','cell_volume'])
fig,ax=plt.subplots(1,1)
ax.step(prof0['the_x'],prof0['the_y'])
ax.step(prof1['the_x'],prof1['the_y'])
axbonk(ax,xlabel=r'$\rho$',ylabel=r'$P(\rho)$',xscale='log',yscale='log')
fig.savefig('/home/dcollins4096/PigPen/PDFs_allcores_n%04d.png'%frame)
fig.savefig('/home/dcollins4096/PigPen/PDFs_allcores_n%04d.pdf'%frame)
