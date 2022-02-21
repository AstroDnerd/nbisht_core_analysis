
from starter2 import *
import matplotlib.image as mpimg

import data_locations as dl
reload(dl)
plt.close('all')
import figure_sublabel
reload(figure_sublabel)
#import three_loopers_1tff as tl
#import three_loopers_mountain_top as TLM

class image_track():
    """here is the doc string"""
    def __init__(self,loop):
        self.this_looper=loop
        self.name=loop.out_prefix
        self.cores_used=[]

    def run(self,core_list=None, external_ax=None, center_image=True):
        dx=1./2048
        nx = 1./dx
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        if external_ax is None:
            fig,ax = plt.subplots(1,1, figsize=(8,8))
        else:
            ax = external_ax
        x_ext = extents()
        y_ext = extents()
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            self.ms = ms
            if ms.nparticles < 10:
                continue
            print('image centroid %s %d '%(self.this_looper.sim_name, core_id))
            self.cores_used.append(core_id)
            self.times = thtr.times

            for LOS in [0]:#,1,2]:
                x = [1,0,0][LOS]
                y = [2,2,1][LOS]
                xlab=r'$%s \rm(code\ length)$'%'xyz'[x]
                ylab=r'$%s \rm(code\ length)$'%'xyz'[y]

                the_pos = np.column_stack( [ms.mean_xc, ms.mean_yc, ms.mean_zc])

                my_x = the_pos[:,x]
                my_y = the_pos[:,y]

                ax.scatter(my_x[0],my_y[0],c='k', s=0.5)
                ax.plot( my_x, my_y, c='k', linewidth=0.3)
                x_ext( my_x)
                y_ext( my_y)


        ax.plot([0,1,1,0,0], [0,0,1,1,0], c=[0.5]*3)

        if center_image:
            dx1 = max( [ np.abs(x_ext.minmax[0]), np.abs( 1-x_ext.minmax[1])])*1.1
            dy1 = max( [ np.abs(y_ext.minmax[0]), np.abs( 1-y_ext.minmax[1])])*1.1
            ddd = max([dx1,dy1])
            ax.set_xlim(-ddd,1+ddd)
            ax.set_ylim(-ddd,1+ddd)
        axbonk(ax,xlabel=xlab,ylabel=ylab)
        figure_sublabel.add_sublabel(ax, self.this_looper)
        ax.set_aspect('equal')
        if external_ax is None:
            output_directory = "./plots_to_sort"
            outname = "%s/centroid_tracks_%s.png"%(output_directory,self.name)
            fig.savefig(outname)
            plt.close(fig)

if 0:
    centroids={}
#for loop in [tl.looper1, tl.looper2, tl.looper3]:
#    name = loop.out_prefix
#    centroids[name] = image_track(loop)
#    centroids[name].run()
    for this_simname in TLM.loops:
        loop = TLM.loops[this_simname]
        name = loop.out_prefix
        centroids[name] = image_track(loop)
        centroids[name].run()

