
from starter2 import *

import convex_hull_tools as CHT
import matplotlib.colors as colors

reload(CHT)
import hair_dryer
reload(hair_dryer)
import stay_close
import three_loopers_tenfour as TL4
sim_list=['u401','u402','u403']
#sim_list=['u402']
import supersets
reload(supersets)
import hair_dryer
reload(hair_dryer)
import colors



class images():
    def __init__(self,this_looper,first_looper):
        self.this_looper=this_looper
        self.first_looper=first_looper
        self.cores_used=[]
        self.name = self.this_looper.sim_name

    def run(self, output_prefix='NAME',core_list=None, frames=None,
            verbose=False, los_list=[-1], external_axis=None):
        print('w3')
        dx=1./2048
        nx = 1./dx
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores
        if hasattr(core_list,'v'):
            core_list=core_list.v #needs to not have unit.
            core_list=core_list.astype('int')
        if frames is None:
            frames = thtr.frames

        thtr.sort_time()

        tsorted = thtr.times/colors.tff
        self.core_list=core_list
        mini_scrubbers = {}
        this_looper=self.this_looper
        if verbose:
            print("MS 1")
        print('w4')
        
        for core_id in core_list:
            ms=  trackage.mini_scrubber(this_looper.tr,core_id, do_velocity=False)
            #ms.make_floats(core_id)

            for nt,frame in [[0,0],[-1,118]]:
                fig,axes=plt.subplots(2,2, figsize=(12,12))
                ax=axes[0][0]; ax1=axes[0][1]
                ax2=axes[1][0];ax3=axes[1][1]
                ax.set_aspect('equal')

                delta = 0.1
                ax.plot([0,1,1,0,0], [0,0,1,1,0], c=[0.5]*3)

                mask=slice(None)

                this_x,this_y,this_z=ms.this_x[mask,nt],ms.this_y[mask,nt], ms.this_z[mask,nt]
                this_p = [this_x,this_y,this_z]

                x=1;y=2

                x_min, x_max, y_min, y_max = [1,0,1,0]

                x_min = min([x_min,this_p[x].min(), -delta])
                x_max = max([x_max,this_p[x].max(), 1+delta])
                y_min = min([y_min,this_p[y].min(), -delta])
                y_max = max([y_max,this_p[y].max(), 1+delta])

                ax.scatter( this_p[x], this_p[y])

                ms2 = trackage.mini_scrubber(self.first_looper.tr, core_id, do_velocity=False)
                first_x,first_y,first_z=ms2.this_x[mask,nt],ms2.this_y[mask,nt], ms2.this_z[mask,nt]
                first_p = [first_x,first_y,first_z]
                ax.scatter( first_p[x], first_p[y])


                other_density = this_looper.tr.c([core_id],'density')[:,nt]
                first_density = self.first_looper.tr.c([core_id],'density')[:,nt]

                other_x = ms.this_x[:,nt]-ms2.mean_x[nt]
                other_y = ms.this_y[:,nt]-ms2.mean_y[nt]
                other_z = ms.this_z[:,nt]-ms2.mean_z[nt]
                other_r = np.sqrt(other_x**2+other_y**2+other_z**2)
                r_ext = extents()
                d_ext = extents()
                ax1.scatter( other_r,other_density, c='r')
                r_ext(other_r); d_ext(other_density)
                ax3.scatter( ms2.r[:,nt], first_density, c='g')
                r_ext(ms2.r[:,nt]); d_ext(first_density)
                axbonk(ax1,xscale='log',yscale='log', xlim=r_ext.minmax,ylim=d_ext.minmax)
                axbonk(ax3,xscale='log',yscale='log', xlim=r_ext.minmax,ylim=d_ext.minmax)



                outname='plots_to_sort/image_%s_%04d.png'%(output_prefix,frame)
                fig.savefig(outname)
                print(outname)


if 0:
    print('w1')
    if 'loopb' not in dir():
        savefile="otherones_b002_temp.h5"
        loopb = looper.core_looper(directory= TL4.loops['u402'].directory,savefile_only_trackage=savefile)
    print('w1')
#hd = hair_dryer.hair_time(loopb)
#color_dict={1:'r'}
#fig,ax=plt.subplots(1,1)
#hd.run(color_dict=color_dict,external_axis=ax)
#fig.savefig('plots_to_sort/hair2.png')

    IM = images(loop, TL4.loops['u402'])
    IM.run(frames=[0,118])
