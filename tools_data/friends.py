from starter2 import *
import annotate_particles_4
reload(annotate_particles_4)


def core_proj_multiple(looper,  field='density', zoom=1,axis_list=[0,1,2], color_dict={},
                    core_list=None,frame_list=None, cores_to_save=None,
                       annotate_grids=True, 
                      annotate_fields=False, annotate_gravity=False, annotate_velocity=False,annotate_velocity_streamlines=False, annotate_lic=False, 
                       annotate_particles=False, annotate_core_ids=True,
                       code_length=True, 
                      tracker_positions=True, shifted_tracker=True, float_positions=False,
                      marker_size=1, verbose=False,
                       zlim=None,
                       cmap='Greys', weight_field=None):
    """
    Plots an collection of cores in a smooth manner.
    FIRST loop over frames,
        Draws boxes around the particles,  
    THEN smooths the path of the edges.
    THEN make all the plots
    """

    tr = looper.tr
    if core_list is None:
        core_list = np.unique(tr.core_ids)
    if cores_to_save is None:
        cores_to_save=core_list
    if frame_list is None:
        frame_list = looper.frame_list
    tracker_index =  [np.where(looper.tr.frames == frame)[0][0] for frame in frame_list]
    times=nar(looper.tr.times[ tracker_index] )
    all_times=looper.tr.times


    #
    #get all the miniscrubbers at once.
    #We should speed this code up.
    #

    if verbose:
        print("Mini scrubbers")
    mini_scrubbers = {}
    for core_id in core_list:
        do_velocity=False
        mini_scrubbers[core_id]=  trackage.mini_scrubber(looper.tr,core_id, do_velocity=do_velocity)
        mini_scrubbers[core_id].make_floats(core_id)

    for frame in frame_list:
        frame_ind = np.where(looper.tr.frames == frame)[0]

        # Check to see if the image was made already,
        # and skips it if it has.
        if len(core_list) == 1:
            suffix = "c%04d"%core_list[0]
        else:
            suffix = 'multi'
        ds = looper.load(frame)

        #
        # main plot loop
        #
        for ax in axis_list:
            proj = ds.proj(field,ax, weight_field=weight_field)
            pw = proj.to_pw(origin='domain')
            
            pw.zoom(zoom)
            if zlim is not None:
                pw.set_zlim(field,zlim[0],zlim[1])

            pw.set_cmap(field,cmap)

            Hcoord=ds.coordinates.x_axis[ax]
            Vcoord=ds.coordinates.y_axis[ax]
            if annotate_lic:
                pw.annotate_line_integral_convolution('magnetic_field_x','magnetic_field_y', lim=(0.5,0.65))
            if annotate_fields:
                #pw.annotate_magnetic_field(plot_args={'color':'b'})
                pw.annotate_streamlines("magnetic_field_x","magnetic_field_y",plot_args={'color':'b'})
            if annotate_velocity:
                pw.annotate_velocity()
            if annotate_velocity_streamlines:
                GH = [YT_velocity_x, YT_velocity_y, YT_velocity_z][Hcoord]
                GV = [YT_velocity_x, YT_velocity_y, YT_velocity_z][Vcoord]
                pw.annotate_streamlines( GH, GV, plot_args={'color':'k'})
            if annotate_gravity:
                GH = [YT_acceleration_x, YT_acceleration_y, YT_acceleration_z][Hcoord]
                GV = [YT_acceleration_x, YT_acceleration_y, YT_acceleration_z][Vcoord]
                pw.annotate_streamlines( GH, GV, plot_args={'color':'g'})


            if annotate_grids:
                pw.annotate_grids()
            if code_length:
                pw.set_axes_unit('code_length')
            if field in ['MagVelAngle']:
                pw.set_cmap('MagVelAngle','hsv')
                pw.set_zlim('MagVelAngle',0,180)
            looper.pw = pw
            if verbose:
                print('save')
            for core_id in core_list:
                ms = mini_scrubbers[core_id]

                positions = np.column_stack([ms.float_x[:,frame_ind], ms.float_y[:,frame_ind], ms.float_z[:,frame_ind]])
                color=color_dict.get(core_id, 'k')
                #color = [1.0,0.0,0.0,0.8]
                #alpha=0.1
                alpha=0.4
                marker_size=1
                if annotate_core_ids:
                    #pw.annotate_sphere(snapshot.R_centroid,Rmax, circle_args={'color':color} ) #R_mag.max())
                    centroid=positions.mean(axis=0)
                    pw.annotate_text(centroid,
                                     "%d"%core_id,text_args={'color':color}, 
                                     inset_box_args={'visible':False},
                                     coord_system='data')
                if annotate_particles:
                    pw.annotate_these_particles4(1.0, col=[color]*positions.shape[0], positions=positions, 
                                                 p_size=marker_size, alpha=alpha)
            for core_id in cores_to_save:
                ms = mini_scrubbers[core_id]
                positions = np.column_stack([ms.float_x[:,frame_ind], ms.float_y[:,frame_ind], ms.float_z[:,frame_ind]])
                centroid=positions.mean(axis=0)
                pw.set_center([ centroid[Hcoord], centroid[Vcoord]])
                print(pw.save( "plots_to_sort/friends_%s_c%04d_n%04d"%(looper.sim_name,core_id,frame)))

sim_list=['u501','u502', 'u503']
import three_loopers_u500 as TL
for sim in sim_list:
    this_looper = TL.loops[sim]
    frame_list = [this_looper.target_frame]
    core_proj_multiple(this_looper,  field='density', zoom=8,axis_list=[0], 
                    core_list=core_list,frame_list=frame_list, 
                       annotate_grids=False, 
                      annotate_fields=False, annotate_gravity=False, annotate_velocity=False,annotate_velocity_streamlines=False, annotate_lic=False, 
                       annotate_particles=False, annotate_core_ids=True,
                       code_length=True, 
                      tracker_positions=True, shifted_tracker=True, float_positions=False,
                      marker_size=1, verbose=False,
                       zlim=None,
                       cmap='Greys', weight_field=None)


