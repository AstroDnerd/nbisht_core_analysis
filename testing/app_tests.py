"""
This displays what all the functions in loop_apps.py do.
"""
from starter2 import *
print(looper)
import xtra_energy
reload(loop_apps)
core_list=[85,86]
frame_list=[0] #[0]# range(0,130,10)
fields=['density']  
#Look for app_test.h5 and read from that.
#If not, make a new looper and write app_test.h5.
save_field = 'all_cores_n0000.h5'
if 'this_looper' not in dir() and os.path.exists(save_field) and True:
    print("WUUUUT")
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    this_looper = looper.core_looper(directory= directory,savefile=save_field)
    this_looper.derived = [xtra_energy.add_force_terms] #for some reason derived quantities aren't saved.

if 'this_looper' not in dir():
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    this_looper = looper.core_looper(directory= directory,
                                     derived=[xtra_energy.add_force_terms],
                                     sim_name = 'u05',
                                     out_prefix = 'plots_to_sort/test',
                                     target_frame = 125,
                                     frame_list = frame_list,
                                     core_list = core_list,
                                     fields_from_grid=['x','y','z']+fields
                                  )
    this_looper.get_target_indices(h5_name='datasets_small/u05_0125_peaklist.h5',
                                     bad_particle_list='datasets_small/bad_particles.h5')
    this_looper.get_tracks()
    #this_looper.save('all_cores_n%04d.h5'%0)

if False: #'this_looper' not in dir():
    directory = '/scratch1/dcollins/Paper42_new_turb/eq45_repeat_eq44/GravPotential'
    this_looper = looper.core_looper(directory= directory,
                                     derived=[],
                                     sim_name = 'eq45',
                                     out_prefix = 'plots_to_sort/test',
                                     target_frame = 1089,
                                     frame_list = [1089],
                                     core_list = list(range(226)),
                                     fields_from_grid=['x','y','z']+fields
                                  )
if 0:   
    this_looper.get_target_indices(h5_name='/scratch1/dcollins/Paper42_new_turb/eq45_repeat_eq44/indices_tile2_000.h5',
                                     bad_particle_list='bad_particles.h5')
if 0:
    this_looper.get_tracks()
    this_looper.save('app_test.h5')

if 0:
    #print the centroid of each core.
    loop_apps.print_centroid(this_looper)

if 0:
    #project the entire domain, with particles form cores in the list plotted
    #All cores together.
    loop_apps.proj_cores(this_looper,axis_list=[0],core_list=[16],field='density')

if 0:
    #this takes a list of particle indices.  
    #mostly for debugging purposes.
    loop_apps.select_particles(this_looper,axis_list=[0],
                               these_particles=this_looper.target_indices[ core_list[0] ])

if 1:
    #Plot each core individually, zoomed in to the core itself.
    loop_apps.core_proj_follow(this_looper,field='density',axis_list=[0],force_log=True)

if 0:
    #Plot each core individually, zoomed in to the core itself.
    this_looper.frame_list=[10]
    loop_apps.slice_raster(this_looper,field='velocity_divergence',axis_list=[0])


if 0:
    """Draws circles around each core."""
    loop_apps.core_circle( this_looper, axis_list=[0])

if 0:
    loop_apps.proj_cores_with_annotations( this_looper, axis_list=[0],color_dict={79:'g',83:'b'})

if 0:
    import testing.early_mask as early_mask
    reload(early_mask)
    for core in [85,86]:
        this_looper.out_prefix='plots_to_sort/u05_c%04d'%core
        early_mask.project_particle_mask(this_looper, axis_list=[2],core_list=[core])

if 0:
    import testing.cic_test as cic
    reload(cic)
    this_looper.core_list = looper.get_all_nonzero()
    print('doing things')
    cic.mask_test4(this_looper, axis_list=[0],core_list=[31])
    #ad=cic.get_deposit_field(this_looper) #,frame=myloop.current_frame, core_list=core_list, mask_stash=mask_stash)
    deposit_tuple=("deposit","deposit_target_particles")
    #this_looper.ad[deposit_tuple]

if 0:
    import testing.cic_test as cic
    reload(cic)
    this_looper.core_list = looper.get_all_nonzero()
    print('doing things')
    cic.mask_test2(this_looper, axis_list=[0])

if 0:
    mask_stash=None
    if 'ad' not in dir():   
        #ds = this_looper.ds
        ds = this_looper.load(frame=this_looper.current_frame,derived=[cic.add_tracer_density])
        cic.add_tracer_density(ds)
        ad = ds.all_data()
    if 1:
        if mask_stash is None:
            mask_stash = np.zeros(ds['NumberOfParticles'], dtype='int32')
        cores = this_looper.target_indices
        cores = this_looper.core_list[0:1]
        all_target_indices = np.concatenate( [this_looper.target_indices[core_id] for core_id in cores])
        mask_stash[:] *= 0 
        my_norm = (8*ds.index.grids[0]['particle_mass'][0]/ds.index.grids[0].dds.prod()).in_units('code_density')
        ad.set_field_parameter('target_indices',all_target_indices)
        ad.set_field_parameter('mask_to_get',mask_stash)
        ad.set_field_parameter('particle_norm',my_norm)
        deposit_tuple=("deposit","deposit_target_particles")
    if 1:
        axis=0
        proj = ds.proj(deposit_tuple,axis,center='c',data_source=ad)
        pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap(deposit_tuple,'gray')
        #outname = '%s_cic_test_%sn%04d'%(looper.out_prefix,snapshot.core_id,snapshot.frame)
        outname = "plots_to_sort/cic_test_%04d"%this_looper.current_frame
        print(outname)
        print( pw.save(outname))
    if 0:
        for g in []: # ds.index.grids[:2]:
            g.set_field_parameter('target_indices',all_target_indices)
            g.set_field_parameter('mask_to_get',mask_stash)
            g.set_field_parameter('particle_norm',my_norm)
            oot=cic.deposit_target_particles_test(None,g)
            #print(g[deposit_tuple].sum())
