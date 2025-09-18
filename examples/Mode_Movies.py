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

if 'mini_scrubbers' not in dir():
    print('mini scrubbers.')
    mini_scrubbers={}
    for sim in ['u501','u502','u503']:
        this_looper=TL.loops[sim]
        mini_scrubbers[sim]={}
        all_cores=np.unique(this_looper.tr.core_ids)
        thtr=this_looper.tr

        for core_id in all_cores:
            print('ms',core_id)
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            ms.particle_pos(core_id)

            mini_scrubbers[sim][core_id]=ms

import find_other_cores
sim_list = ['u501','u502','u503']
sim_list = ['u503']
for sim in sim_list:
    set_looper = TL.loops[sim]
    main_core_list = np.unique(TL.loops[sim].tr.core_ids)
    main_core_list.sort()
    main_core_list=main_core_list[::-1]
    #main_core_list = [62]
    other_cores= find_other_cores.get_other_cores( set_looper, main_core_list, mini_scrubbers[sim])

    if 1:
        movie_mask = movie_frames.quantized_mask(set_looper)
        frame_list=set_looper.tr.frames[movie_mask.flatten()]
        #frame_list = frame_list[ frame_list > 71]


    #camera = camera_path.camera_1( set_looper, 'tight_float')
    #camera = camera_path.camera_1( set_looper, 'domain')
    #camera = camera_path.camera_1( set_looper, 'smooth_zoom')
    camera = camera_path.camera_1( set_looper, 'smooth_zoom_2')
    


    #frame_list=[0, set_looper.target_frame] #frame_list[::10]
    set_looper.plot_directory = "./MOVIE_PLOTS"
    #set_looper.plot_directory = "./plots_to_sort"
    cores_yet_to_do=[]
    for main_core in main_core_list:
        print("MAIN",sim,main_core)

        color_dict = colors.make_core_cmap(other_cores[main_core], cmap = 'autumn', seed = -1)
        color_dict[main_core]=[0.0,1.0,0.0]
        core_list = sorted(other_cores[main_core]+[main_core])

        

        set_looper.ds_list={}
        derived=[]
        #derived=[xtra_energy.add_b_over_rho]
        derived=[xtra_energy.add_energies]
        derived=[ xtra_energy.add_energies, xtra_energy.add_gravity]
        set_looper.derived=derived
        mono=core_proj.core_proj_multiple(set_looper,camera=camera,axis_list=[0,1,2], 
                                          main_core=main_core,
                                          clobber=False,
                                          cmap='Greys',
                                     color_dict=color_dict,
                                     frame_list = frame_list,
                                     field=YT_density,
                                          velocity=False,
                                          mean_velocity=False,
                                     core_list=core_list,
                                     slab=False, zoom=True, only_sphere=True, center_on_sphere=True, 
                                     annotate=True,  plot_particles=True,
                                     grids=False, float_positions=True, monotonic=False, 
                                          verbose=True,
                                          path_only=False, 
                                         plot_y_tracks=False, plot_points=False,
                                         derived=derived, zlim=[1e-3,1e3])
                                   #zoom=False, grids=False, particles=False, moving_center=True) #, frame_list=[1])
#loop_apps.core_proj_multiple(this_looper,axis_list=[0], field='density', color_dict=color_dict, particles=True, fields=False, core_list )#, frame_list=[31])
