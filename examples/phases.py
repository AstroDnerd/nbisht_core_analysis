from starter2 import *
import xtra_energy
import data_locations as dl
import loop_tools
reload(dl)

reload(looper)
reload(loop_apps)
reload(loop_tools)
import three_loopers_u500 as TL5
import colors
import core_proj
import movie_frames
reload(movie_frames)
import camera_path
reload(camera_path)
reload(core_proj)
sim_list = ['u501','u502','u503']
sim_list = ['u502']




for this_simname in sim_list:
    set_looper = TL5.loops[this_simname]

    all_cores=  np.unique(set_looper.core_ids)


#frame_list=set_looper.tr.frames #[::10]
    if 1:
        movie_mask = movie_frames.quantized_mask(set_looper)
        frame_list=set_looper.tr.frames[movie_mask]
        #frame_list = frame_list[ frame_list > 71]


    camera = camera_path.camera_1( set_looper, 'tight_float')
#camera = camera_path.camera_1( set_looper, 'domain')
#camera = camera_path.camera_1( set_looper, 'smooth_zoom')
    camera = camera_path.camera_1( set_looper, 'smooth_zoom_2')

#frame_list=frame_list[::10]
#frame_list=[frame_list[0],frame_list[-1]]
#color_dict = colors.make_core_cmap(core_list, cmap = 'tab20c', seed = -1)
#frame_list = [0,50,111]
    #frame_list=[0,1,2]
    set_looper.plot_directory = "./plots_to_sort"
    frame_list=frame_list[::5]
    import phase_movie
    reload(phase_movie)
    for core_id in [74]:
        core_list=[core_id]
        color_dict={core_id:'k'}
        set_looper.ds_list={}
        derived=[]
        #derived=[xtra_energy.add_b_over_rho]
        derived=[xtra_energy.add_energies]
        derived=[ xtra_energy.add_energies, xtra_energy.add_gravity]
        set_looper.derived=derived
        mono=phase_movie.phase_movie(set_looper,camera=camera,
                                          fields=[YT_grav_energy_2,YT_kinetic_energy],
                                          frame_list=frame_list,core_list=core_list,only_sphere=True)
