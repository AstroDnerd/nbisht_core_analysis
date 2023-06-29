
from starter2 import *

import core_proj_three
reload(core_proj_three)

import camera_path
import track_loader as TL
sim_list=['u501','u502','u503']
sim_list=['u502']
sim_list=['t002']
TL.load_tracks(sim_list)
import movie_frames
if 1:
    for sim in sim_list:
        loop = TL.loops[sim]
        camera = camera_path.camera_1( loop, 'smooth_zoom_2')
        frames = nar(loop.frame_list)
        #mask = movie_frames.quantized_mask(loop)
        #mask[-1]=True
        #frames = frames[mask]
        if 0:
            every_ten = np.zeros_like(frames,dtype='bool')
            every_ten[::10]=True
            every_ten[-1]=True
            frames = frames[every_ten]

        core_proj_three.core_proj_multiple(loop,
                                          # axis_list=[0],
                                           core_list=[0],frame_list=frames,camera=camera, main_core=0,
                                           only_sphere=True,
                                          monotonic=False)
