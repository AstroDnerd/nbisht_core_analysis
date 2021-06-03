"""
This displays what all the functions in loop_apps.py do.
"""
from starter2 import *
import xtra_energy
reload(loop_apps)
#core_list=[0,1, 10, 27,79, 16, 70, 44]
core_list = looper.get_all_nonzero() 
frame_list=[125] #range(0,130,10)
fields=['density']  


def proj_cores_annotate_zoom(this_looper, axis_list=[0,1,2],core_list=None, center_cores=[],field='density',color='r',
                            zoom_level=4,cb_label=None, frame_list=None,plot_dir="./plots_to_sort", annotate_particles=False):
    """Full projections of the data, with core particles marked."""
    if frame_list is None:
        frame_list = [this_looper.target_frame]
    if core_list is None:
        core_list = list(this_looper.target_indices.keys())
    if core_list == "short":
        core_list = list(this_looper.target_indices.keys())[:5]
    if center_cores == None:
        center_cores = core_list
    if center_cores == "one":
        center_cores = core_list[0:1]
    do_peak_proj=False
    if center_cores == [] or center_cores == 'center':
        #make sure it saves
        do_peak_proj=True
    do_center_proj = False
    if center_cores == 'center':
        center_cores = []
        do_center_proj = True

    for frame in frame_list:
        print("%s %d"%(this_looper.out_prefix,frame))
        for snapshot in this_looper.snaps[frame].values():
            if snapshot.R_centroid is None:
                snapshot.get_all_properties()
        for axis in axis_list:
            ds = this_looper.load(frame)
            center = this_looper.ds.arr(nar([0.5]*3),'code_length')
            #this_looper.proj = yt.ProjectionPlot(this_looper.ds, axis=axis, fields=[field])#, center=center)
            proj = ds.proj(field,axis,center=center)
            this_looper.proj=proj.to_pw(width=(1.0,'code_length'),center=center,origin='domain')
            #this_looper.proj.set_cmap(field,'Greys')
            this_looper.proj.set_cmap(field,'Greys')
            this_looper.proj.zoom(zoom_level)
            this_looper.proj.set_axes_unit('code_length')
            for nc,core_number in enumerate(core_list):
                snapshot = this_looper.snaps[frame][core_number]
                center = ds.arr(snapshot.R_centroid,'code_length')
                this_looper.proj.annotate_text(loop_apps.cbump(center),
                                 ".",text_args={'color':color}, 
                                 inset_box_args={'visible':False},
                                 coord_system='data')
                this_looper.proj.annotate_text(loop_apps.cbump(center),
                                 "%d"%snapshot.core_id,text_args={'color':color}, 
                                 inset_box_args={'visible':False},
                                 coord_system='data')
                if annotate_particles:
                    this_looper.proj.annotate_select_particles4(1.0, col='r', indices=this_looper.target_indices[core_number])
            if do_center_proj or do_peak_proj:
                outname = '%s/%s_core_zoom_annotate_n%04d'%(plot_dir,this_looper.out_prefix,frame)
                if do_center_proj:
                    center = ds.arr([0.5]*3,'code_length')
                else:
                    maxval, center = ds.find_max('density')
                plot_x = [1,2,0][axis]
                plot_y = [2,0,1][axis]
                this_center = center[plot_x],center[plot_y]
                this_looper.proj.set_center(this_center)
                print( this_looper.proj.save(outname))

            for core_number in center_cores:
                snapshot = this_looper.snaps[frame][core_number]
                center = ds.arr(snapshot.R_centroid,'code_length')
                plot_x = [1,2,0][axis]
                plot_y = [2,0,1][axis]
                this_center = center[plot_x],center[plot_y]
                outname = '%s/%s_core_zoom_annotate_c%04d_n%04d'%(plot_dir,this_looper.out_prefix,core_number,frame)
                if cb_label is not None:
                    this_looper.proj.set_colorbar_label(field,cb_label)
                this_looper.proj.set_center(this_center)
                print( this_looper.proj.save(outname))

import three_loopers_1tff as tl
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

