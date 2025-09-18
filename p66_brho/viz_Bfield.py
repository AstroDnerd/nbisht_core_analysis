
from starter2 import *
import davetools
reload(davetools)
import p49_fields
reload(p49_fields)
import math
import pcolormesh_helper as pch

from mpl_toolkits.mplot3d import Axes3D
from yt.visualization.api import Streamlines
# --- --- --- --- --- --- ---

class visualizeB(): 
    def __init__(self,the_loop):
        self.this_looper = the_loop

    def qtyRun(self,core_list=None): 
        print('inside qtyRun')
        thtr = self.this_looper.tr

        # CORES
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        # THE FRAMES
        #the_frame = thtr.frames[-1:] 
        # EVERY TEN FRAMES - make a loop now!
        for nf,frame in enumerate(thtr.frames): 

            # CORE-LOOP
            for nc,core_id in enumerate(core_list):
                ds = self.this_looper.load(frame) 
                ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
                ms.particle_pos(core_id)

                the_center = ms.mean_center[:,-1]  #the three coords for the last frame 
                the_radius = 1/128  
                diameter = the_radius * 2
                the_sphere = ds.sphere(the_center, the_radius) 
                
                if 0:
                    crosssect = yt.SlicePlot(ds, "z", ("gas", "density"), center="c", width=diameter)
                    crosssect.annotate_streamlines(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"))
                    crosssect.save()

                # plottingn streamlines
                N = 100 #number of streamlines
                scale = diameter
                pos_dx = np.random.random((N, 3)) * scale - scale / 2.0
                pos = the_center + pos_dx  #needs pos_dx

                print('before streamlines')
                b_streamlines = Streamlines(ds, pos, ("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"), ("gas", "magnetic_field_z"),\
                                            get_magnitude = True) 
                print('after streamlines')
                b_streamlines.integrate_through_volume()

                fig = plt.figure()
                ax = Axes3D(fig, auto_add_to_figure=False)
                fig.add_axes(ax)

                for stream in b_streamlines.streamlines:
                    stream = stream[np.all(stream != 0.0, axis=1)]
                    #pdb.set_trace() 
                    ax.plot3D(stream[:, 0], stream[:, 1], stream[:, 2], alpha=0.4)

                # Save the plot to disk.
                plt.savefig("streamlines_%d_%d_602.png"%(core_id,frame))


# MAIN
import three_loopers_six as TL6
if 'viz' not in dir():
    viz=True
if 'viz2' not in dir() or viz:
    viz2=visualizeB(TL6.loops['u602'])

if 1:
    all_cores = np.unique(viz2.this_looper.tr.core_ids)
    #core_list = all_cores
    core_list = all_cores[:1]  #DEBUG
    #core_list = core_list1[:3] 
    print('before running')
    viz2.qtyRun(core_list=core_list) 





