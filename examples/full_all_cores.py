from starter2 import *
import xtra_energy
import data_locations as dl
import loop_tools
reload(dl)
import camera_path
reload(camera_path)

reload(looper)
reload(loop_apps)
reload(loop_tools)
import colors
import core_proj
reload(core_proj)
if 'this_simname' not in dir():
    this_simname = 'u503'
import three_loopers_u500 as TL5
set_looper = TL5.loops[this_simname]

all_cores=  np.unique(set_looper.core_ids)


if 1:
    import movie_frames
    reload(movie_frames)
    movie_mask = movie_frames.quantized_mask(set_looper)
    frame_list=set_looper.tr.frames[movie_mask]
    #frame_list = frame_list[ frame_list > 71]


#camera = camera_path.camera_1( set_looper, 'tight_float')
camera = camera_path.camera_1( set_looper, 'domain')
#camera = camera_path.camera_1( set_looper, 'smooth_zoom')
#camera = camera_path.camera_1( set_looper, 'smooth_zoom_2')

#frame_list=frame_list[::10]
#frame_list=[frame_list[0],frame_list[-1]]
#color_dict = colors.make_core_cmap(core_list, cmap = 'tab20c', seed = -1)
#frame_list=[111]
if 1:
    core_list=None #[16]
    color_dict={}
    set_looper.ds_list={}
    derived=[]
    #derived=[xtra_energy.add_b_over_rho]
    derived=[xtra_energy.add_energies]
    derived=[ xtra_energy.add_energies, xtra_energy.add_gravity]
    set_looper.derived=derived
    set_looper.out_prefix = 'full_all_centers_%s'%set_looper.sim_name
    mono=core_proj.core_proj_multiple(set_looper,camera=camera,axis_list=[0], 
                                      #cmap='seismic',
                                      #field=YT_ge_ke,
                                      # zlim=[1e-3,1e3],
                                      cmap='Greys',
                                      field=YT_density,
                                      color_dict=color_dict,
                                      frame_list = frame_list,
                                      velocity=False,
                                      mean_velocity=False,
                                      core_list=core_list,
                                      slab=False, zoom=True, only_sphere=False, center_on_sphere=False, 
                                      annotate=False,  
                                      plot_particles=False,
                                      grids=False, float_positions=True, monotonic=False, 
                                      verbose=True,
                                      path_only=False, 
                                      plot_y_tracks=False, plot_points=False,
                                      derived=derived)
                               #zoom=False, grids=False, particles=False, moving_center=True) #, frame_list=[1])
#loop_apps.core_proj_multiple(this_looper,axis_list=[0], field='density', color_dict=color_dict, particles=True, fields=False, core_list )#, frame_list=[31])
