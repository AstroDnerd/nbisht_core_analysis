
from starter2 import *

reload(loop_apps)

#import three_loopers
make_new_looper=False
if 'this_simname' not in dir():
    this_simname = 'u11'
if make_new_looper:
    all_nonzero = looper.get_all_nonzero(dl.n_particles[this_simname])
    core_list = all_nonzero.astype('int')[::-1]
    target_frame = dl.target_frames[this_simname]
    frame_list = [0,80]# [0]+list(range(10,target_frame,10))+[target_frame]
    fields = ['x','y','z','density']
    fields += ['vorticity_magnitude']
    derived=[]
    print(core_list)
    new_looper = looper.core_looper(directory= dl.sims[this_simname],
                                     sim_name = this_simname,
                                     out_prefix = this_simname,
                                     target_frame = dl.target_frames[this_simname],
                                     frame_list = frame_list,
                                     core_list =  core_list,
                                     fields_from_grid=fields,
                                     derived = derived
                                  )
    print(core_list)
    new_looper.get_target_indices(h5_name=dl.peak_list[this_simname],
                                     bad_particle_list=dl.bad_particles.get(this_simname,None))
    new_looper.get_tracks()
    new_looper.save('temp.h5')

for frame in [0,80]:
#    loop_apps.phase_with_preimage(looper2,frame,['density','magnetic_field_strength','cell_volume'],weight_field=None,
#                                 xlim=(1e-3,1e7),ylim=[1e-1,1e4],zlim=[1e-10,1e-1])
    loop_apps.phase_with_preimage(new_looper,frame,['density','vorticity_magnitude','cell_volume'],weight_field=None,
                                  xlim=(1e-3,1e7),ylim=[1e-1,5e4])#,zlim=[1e-1,1e4])
