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


def proj_cores_annotate_zoom(this_looper, axis_list=[0,1,2],core_list=None, center_cores=None,field='density',color='r',
                            zoom_level=4,cb_label=None, frame_list=None,plot_dir="./plots_to_sort", annotate_particles=False):
    """Full projections of the data, with core particles marked."""
    if frame_list is None:
        frame_list = this_looper.frame_list
    if core_list is None:
        core_list = this_looper.core_list
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
            position_dict={}
            center_dict={}
            for nc,core_number in enumerate(core_list):
                ms = trackage.mini_scrubber(this_looper.tr,core_number, do_velocity=False)
                frame_ind = np.where(this_looper.tr.frames == frame)[0][0]
                shifted_tracker=True
                if shifted_tracker:
                    this_x=ds.arr(ms.this_x[:,frame_ind],"code_length")
                    this_y=ds.arr(ms.this_y[:,frame_ind],"code_length")
                    this_z=ds.arr(ms.this_z[:,frame_ind],"code_length")
                else:
                    this_x=ds.arr(ms.raw_x[:,frame_ind],"code_length")
                    this_y=ds.arr(ms.raw_y[:,frame_ind],"code_length")
                    this_z=ds.arr(ms.raw_z[:,frame_ind],"code_length")

                positions = np.column_stack([this_x,this_y,this_z])
                position_dict[core_number]=positions
                center = ds.arr(ms.mean_center_density[:,frame_ind],'code_length')
                center_dict[core_number]=center
                this_looper.proj.annotate_text(center,
                                 ".",text_args={'color':color}, 
                                 inset_box_args={'visible':False},
                                 coord_system='data')
                this_looper.proj.annotate_text(center,
                                 "%d"%core_number,text_args={'color':color}, 
                                 inset_box_args={'visible':False},
                                 coord_system='data')
                if annotate_particles:
                    #this_looper.proj.annotate_select_particles4(1.0, col='r', indices=this_looper.target_indices[core_number])
                    this_looper.proj.annotate_these_particles2(1.0, col=[color]*positions.shape[0], positions=positions)
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
                center = center_dict[core_number]
                print("CENTER",center)
                plot_x = [1,2,0][axis]
                plot_y = [2,0,1][axis]
                this_center = ds.arr([center[plot_x],center[plot_y]], 'code_length')
                outname = '%s/%s_core_zoom_annotate_c%04d_n%04d'%(plot_dir,this_looper.out_prefix,core_number,frame)
                if cb_label is not None:
                    this_looper.proj.set_colorbar_label(field,cb_label)
                print(this_center)
                this_looper.proj.set_center(this_center)
                print( this_looper.proj.save(outname))

