from starter2 import *
import xtra_energy
import data_locations as dl
import loop_tools
reload(dl)

reload(looper)
reload(loop_apps)
reload(loop_tools)

this_simname = 'u301'

import three_loopers_mountain_top as TLM
this_looper = TLM.loops[this_simname]


if 'this_looper' not in dir():
    output_base = "%s_cores"%this_simname
#Cores that overlap more than 90% with u203 core 84
    core_list =  [323]
#core_list = [84, 112]
    frame_list = list(range(0,100,10))+list(range(100,dl.target_frames[this_simname]+1))
    fields = ['x','y','z','density']
    derived = []
    #for making
    this_looper = looper.core_looper(directory= dl.sims[this_simname],
                                     sim_name = this_simname,
                                     out_prefix = this_simname,
                                     target_frame = dl.target_frames[this_simname],
                                     frame_list = frame_list,
                                     core_list =  core_list,
                                     fields_from_grid=fields,
                                     derived = derived,
                                     do_shift=False
                                  )
    this_looper.plot_directory = "/home/dccollins/PigPen"
    this_looper.get_target_indices(h5_name=dl.peak_list[this_simname])
    #                                 bad_particle_list=dl.bad_particles.get(this_simname,None))
    this_looper.get_tracks()
    #this_looper.make_snapshots()

import colors
#color_dict = colors.make_core_cmap( this_looper.core_list)
color_dict={323:'r'}
#loop_apps.core_proj_follow(this_looper,axis_list=[0], field='PotentialField',
#                           zoom=False, grids=False, particles=False, moving_center=True) #, frame_list=[1])
#slab = {'zmin':zmin, 'zmax':zmax}
#loop_apps.core_proj_follow(this_looper,axis_list=[2], field='PotentialField', slab=slab, zoom=False, only_sphere=False, center_on_sphere=False)#, frame_list=[31])
##                           #zoom=False, grids=False, particles=False, moving_center=True) #, frame_list=[1])
loop_apps.core_proj_multiple(this_looper,axis_list=[0], field='density', 
                             core_list = [323],
                             #frame_list = [100],
                             zoom=True,
                             only_sphere=False,
                             marker_size=7,
                             color_dict=color_dict, particles=True, fields=False,
                             float_positions=True)
