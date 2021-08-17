
import three_loopers_1tff as tl
import projections
reload(projections)
proj_cores_annotate_zoom = projections.proj_cores_annotate_zoom
if 0:
    #full projections
    for this_looper in [tl.looper1,tl.looper2,tl.looper3]:
        proj_cores_annotate_zoom(this_looper,axis_list=[0],core_list=[],center_cores=[],cb_label='density',
                                plot_dir="./plots_to_sort", zoom_level=1)
if 0:
    #Core 8 and its neighbors
    core_list = [7,8,9,12,252,251,253,10,11,248,248,247,250,133,135]
    center_cores = [8]
    for this_looper in [tl.looper1]:
        proj_cores_annotate_zoom(this_looper,axis_list=[0],core_list=core_list,center_cores=center_cores,cb_label='density',
                                plot_dir="./plots_to_sort", zoom_level=4)
if 1:
    #Core 8 and its neighbors: pre-image, with particles
    core_list = [7,8,9,12,252,251,253,10,11,248,248,247,250,133,135]
    center_cores = [8]
    for this_looper in [tl.looper1]:
        proj_cores_annotate_zoom(this_looper,axis_list=[0],core_list=core_list,center_cores='center',cb_label='density',
                                plot_dir="./plots_to_sort", zoom_level=1, annotate_particles=True, frame_list=[1])

if 0:

    for this_looper in [tl.looper1,tl.looper2,tl.looper3]:
        proj_cores_annotate_zoom(this_looper,axis_list=[0],core_list=None,center_cores=None,cb_label='density',
                                plot_dir="./plots_to_sort_2")

if 0:
    proj_cores_annotate_zoom(this_looper,axis_list=[0],core_list=core_list,center_cores=[34],cb_label='density',
                             zoom_level=8)

if 0:
    #project the entire domain, with particles form cores in the list plotted
    #All cores together.
    if 0:
        this_looper.get_target_indices(h5_name='datasets_small/u05_0125_peaklist.h5',
                                         bad_particle_list='bad_particles.h5')
        this_looper.get_tracks()
    loop_apps.proj_cores2(this_looper,axis_list=[0],core_list=core_list,field='density')

if 0:
    import  datasets_small.u05_core_speeds as sped
    color_dict = {}
    this_looper.frame_list=list(range(0,125,10))+[125]
    for n, regime in enumerate([ [sped.fast_cores,'r'], [sped.ok_cores,'g'], 
                                [sped.slow_cores,'b'], [sped.small_cores, 'c']]):
        color = regime[1]
        this_dict = dict( zip(regime[0], [color]*len(regime[0])))
        color_dict.update(this_dict)

    this_looper.out_prefix='plots_to_sort/proj_with_regimes'
    loop_apps.proj_with_species(this_looper,axis_list=[0],core_list=core_list,field='density', color_dict=color_dict)

