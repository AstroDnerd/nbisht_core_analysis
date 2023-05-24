from starter2 import *

save_field = '../Datasets/all_primitives/all_primitives_c0016.h5'
file_list = glob.glob('../Datasets/all_primitives/all_primitives_c*.h5')
if 'tll' not in dir():# and os.path.exists(save_field):
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    tll = looper.core_looper(directory= directory)#,savefile=save_field)
    for nfile,fname in enumerate(file_list):
        tll.load_loop(fname)
        print(tll.core_list)
        print( "Reading file %d of %d"%(nfile,len(file_list)))
    thtr = tll.tr
    thtr.sort_time()
    all_cores = np.unique(thtr.core_ids)

if 'tll' not in dir() and False:
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    frame_list=list(range(50,60))
    core_list=[16]
    fields=['density']
    tll = looper.core_looper(directory= directory,
                                     derived=[],#xtra_energy.add_force_terms],
                                     sim_name = 'u05',
                                     out_prefix = 'plots_to_sort/follow',
                                     target_frame = 125,
                                     frame_list = frame_list,
                                     core_list = core_list,
                                     fields_from_grid=['x','y','z']+fields
                                  )
    tll.get_target_indices(h5_name='datasets_small/u05_0125_peaklist.h5',
                                     bad_particle_list='datasets_small/bad_particles.h5')
    tll.get_tracks()
