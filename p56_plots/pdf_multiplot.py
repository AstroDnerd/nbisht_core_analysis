from starter2 import *

color={'u05':'r','u10':'g','u11':'b'}
out_dir='plots_to_sort'
these_simnames=['u05','u10','u11']

class prof_from_disk():
    """there is probably a different reader"""
    def __init__(self,profname):
        self.fname = profname
        self.d = {}
        self.read(self.fname)

    def read(self,fname):
        fptr = h5py.File(fname,'r')
        for a in fptr['data']:
            self.d[a]=fptr['data'][a][()]

        fptr.close()



prof_all = {}
prof_mask={}
field = 'potential'
field='density'
for this_simname in these_simnames:
    prof_all[this_simname]  = prof_from_disk("%s_%s_pdf_all_n0000.h5"%(this_simname,field))
    prof_mask[this_simname] = prof_from_disk("%s_%s_pdf_mask_n0000.h5"%(this_simname,field))


fig,ax=plt.subplots(1,1)
for this_simname in these_simnames:
    c=color[this_simname]
    this_prof = prof_all[this_simname]
    ax.plot( this_prof.d['x'], this_prof.d['cell_volume'],label=this_simname,c=c)
    this_prof = prof_mask[this_simname]
    ax.plot( this_prof.d['x'], this_prof.d['target_particle_volume'],label=this_simname,linestyle='--',c=c)

if field == 'potential':
    axbonk(ax,xlabel=r'$\phi$',ylabel=r'$V(\phi)$',yscale='log')
elif field == 'density':
    axbonk(ax,xlabel=r'$\rho$',ylabel=r'$V(\rho)$',yscale='log',xscale='log')

fig.savefig('%s/%s_%s_multiplot.pdf'%(out_dir,'uXX',field))
plt.close(fig)
