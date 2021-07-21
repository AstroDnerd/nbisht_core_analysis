
import matplotlib.colors as colors
from starter2 import *
import matplotlib.image as mpimg

from scipy.optimize import curve_fit
import data_locations as dl
from collections import defaultdict
reload(dl)
plt.close('all')



class rho_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=nar([])
        self.slopes=defaultdict(list)
        self.times=[]

    def run(self,do_all_plots=True,core_list=None,plot_each_frame=False):


        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        self.times=thtr.times[ thtr.times > 0] 
        fig, ax = plt.subplots(1,1)
        for core_id in core_list:
            density = thtr.c([core_id],'density')
            ax.scatter(density[:,0], density[:,-1])
        fig.savefig('plots_to_sort/%s_density_density.pdf'%self.this_looper.out_prefix)

import three_loopers as TL
if 'giant_looper' not in dir() or clobber:
    giant_looper={}
    for loo in [TL.looper1]:
        ds_2 = loo.load( loo.target_frame)
        prefix = loo.out_prefix
        ds_0 = loo.load( 0 )
        ad_2 = ds_2.all_data()
        ad_0 = ds_0.all_data()
        all_particles = ad_0['particle_index'][::100]
        fields=['density']
        giant_looper[prefix] = looper.core_looper(directory=loo.directory,
                                                  sim_name = prefix, out_prefix=prefix,
                                                  target_frame=loo.target_frame,
                                                  frame_list=[0,loo.target_frame],
                                                  core_list=[0],
                                                  fields_from_grid=fields)
        giant_looper[prefix].target_indices[0]=all_particles
        giant_looper[prefix].get_tracks()
for prefix in giant_looper:
    loo = giant_looper[prefix]
    thtr = giant_looper[prefix].tr
    d0 = thtr.track_dict['density'][:,0]
    d1 = thtr.track_dict['density'][:,1]

    fig,ax=plt.subplots(1,1)
    h,xedge,yedge=np.histogram2d(np.log10(d0), np.log10(d1))
    norm = colors.LogNorm(vmin=1,vmax=h.max())
    plot=ax.imshow(h,norm=norm, interpolation='nearest',origin='lower')
    cb1 = fig.colorbar(plot)
    fig.savefig('plots_to_sort/%s_density_density.pdf'%prefix)






#rt1 = rho_tool(TL.looper1)
#rt1.run()
