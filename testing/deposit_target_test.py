from starter2 import *


import testing.early_mask as em
reload(em)

if 'loop_dict' not in dir():
    if 0:
        #take 1, not enough particles.
        import three_loopers_1tff as tl
        loop_dict = {'u201':tl.looper1, 'u202':tl.looper2, 'u203':tl.looper3}
    elif 0:
        #take 1.5, a slightly different core set
        loop_dict = {}
        core_list={}
        core_set={'u05':'/data/cb1/Projects/P19_CoreSimulations/CoreSets/all_cores_n0000.h5',
                  'u10':'/data/cb1/Projects/P19_CoreSimulations/CoreSets/u10_primitives_cXXXX_n0000.h5',
                  'u11':'/data/cb1/Projects/P19_CoreSimulations/CoreSets/u11_primitives_cXXXX_n0000.h5'}
        for this_simname in ['u05','u10','u11']:
            directory = dl.sims[this_simname]
            save_field = core_set[this_simname]
            print("LOAD", this_simname)
            loop_dict[this_simname] = looper.core_looper(directory= directory,savefile=save_field)
            core_list[this_simname] =looper.get_all_nonzero(dl.n_particles[this_simname])
    else:
        import three_loopers_mountain_top as TLM
        reload(TLM)
        loop_dict = TLM.loops

this_looper = loop_dict['u301']
frame = 0
if 0:
    all_target_indices = np.concatenate( [this_looper.target_indices[core_id] for core_id in this_looper.core_list])
if 1:
    core_list=[10,32]
    core_list = this_looper.core_list
    all_target_indices = np.concatenate( [this_looper.target_indices[core_id] for core_id in core_list])
    all_target_indices = all_target_indices.astype('int64')

if 1:
    ds = this_looper.load(frame=frame,derived=[em.add_tracer_density])
    em.add_tracer_density(ds)
    ad = ds.all_data() #ds.region([0.5]*3,[0.4]*3,[0.6]*3)
    ad.set_field_parameter('target_indices',all_target_indices)
    ad.set_field_parameter('mask_to_get',np.zeros_like(all_target_indices,dtype='int32'))
    deposit_tuple = ("deposit","target_particle_volume")

    mask = (ad[deposit_tuple] > 0 )
    density_m = ad['density'][mask]
    density_t = this_looper.tr.c(core_list,'density')[:,0]
    v_m = np.sort( ad['velocity_magnitude'][mask] )
    v_t = np.sort( this_looper.tr.c(core_list,'velocity_magnitude')[:,0] )
    #v_t_2 = np.sort( this_looper.tr.track_dict['velocity_magnitude'][:,0])
