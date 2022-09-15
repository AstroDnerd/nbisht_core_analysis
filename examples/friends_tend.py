from starter2 import *
import xtra_energy
import data_locations as dl
import loop_tools
import movie_frames
reload(movie_frames)
import camera_path
reload(camera_path)
reload(dl)

reload(looper)
reload(loop_apps)
reload(loop_tools)
import three_loopers_u500 as TL
import core_proj
reload(core_proj)

import find_other_cores
sim_list = ['u501','u502','u503']
sim_list = ['u502']
for sim in sim_list:
    set_looper = TL.loops[sim]
    frame_list=[set_looper.target_frame]
    core_list=sorted(np.unique(set_looper.tr.core_ids))
    cores_to_center = [193,74,112,367]


    #camera = camera_path.camera_1( set_looper, 'tight_float')
    camera = camera_path.camera_1( set_looper, 'domain')
    #camera = camera_path.camera_1( set_looper, 'smooth_zoom')
    #camera = camera_path.camera_1( set_looper, 'smooth_zoom_2')
    
    #frame_list=[0, set_looper.target_frame] #frame_list[::10]
    #set_looper.plot_directory = "./MOVIE_PLOTS"
    set_looper.plot_directory = "./plots_to_sort"



    reload(core_proj)
    set_looper.ds_list={}
    derived=[]
    #derived=[xtra_energy.add_b_over_rho]
    derived=[xtra_energy.add_energies]
    derived=[ xtra_energy.add_energies, xtra_energy.add_gravity]
    set_looper.derived=derived
    mono=core_proj.core_proj_multiple(set_looper,camera=camera,axis_list=[0,1,2], 
                                      main_core=main_core,
                                      clobber=True,
                                      cmap='Greys',
                                 frame_list = frame_list,
                                 field=YT_density,
                                      velocity=False,
                                      mean_velocity=False,
                                 core_list=core_list,
                                 slab=False, zoom=None, only_sphere=False, center_on_sphere=True, 
                                 annotate=True,  plot_particles=False,
                                 grids=False, float_positions=True, monotonic=False, 
                                      verbose=True, no_margin=True,
                                      path_only=False, 
                                     plot_y_tracks=False, plot_points=False, cores_to_center=cores_to_center)
                               #zoom=False, grids=False, particles=False, moving_center=True) #, frame_list=[1])
#loop_apps.core_proj_multiple(this_looper,axis_list=[0], field='density', color_dict=color_dict, particles=True, fields=False, core_list )#, frame_list=[31])
