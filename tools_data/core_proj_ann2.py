from starter2 import *
import annotate_particles_4
reload(annotate_particles_4)


class MonotoneColorbar():
    def __init__(self):
        self.prev=None
    def __call__(self,val, verbose=False):
        new_val = val+0
        if self.prev is not None:
            new_val[0] = min( [self.prev[0], min(new_val)])
            new_val[1] = max( [self.prev[1], max(new_val)])
            print(val, new_val)
        self.prev = new_val
        return new_val


def core_proj_multiple(looper, camera=None, field='density', axis_list=[0,1,2], color_dict={},force_log=None,linthresh=100,
                       main_core=None, mini_scrubbers=None,
                    core_list=None,frame_list=None, clobber=True, weight_field=None,
                       only_sphere=True, 
                       zoom=True, moving_center=False, slab=None,
                       grids=True, plot_particles=True, annotate=False, 
                      fields=False, velocity=False, mean_velocity=False, lic=False, 
                       code_length=True, 
                      tracker_positions=True, shifted_tracker=True, monotonic=False, float_positions=False,
                      marker_size=1, path_only=False, verbose=False, plot_y_tracks=True, plot_points=False,
                       zlim=None, no_margin=False, cores_to_center=[],
                      derived=[], cmap='Greys'):
    """
    Plots an collection of cores in a smooth manner.
    FIRST loop over frames,
        Draws boxes around the particles,  
    THEN smooths the path of the edges.
    THEN make all the plots
    """

    if core_list is None:
        core_list = looper.core_list
    if frame_list is None:
        frame_list = looper.frame_list
    all_png = glob.glob("%s/*png"%looper.plot_directory)
    tr = looper.tr
    if monotonic is True:
        #monotonic  ={'zlim':MonotoneEnforcer()}#, 'left':MonotoneEnforcer2(nvalues=6),'right':MonotoneEnforcer2(nvalues=6)}
        monotonic  ={'zlim':MonotoneColorbar()}#, 'left':MonotoneEnforcer2(nvalues=6),'right':MonotoneEnforcer2(nvalues=6)}
        #monotonic  ={'zlim':MonotoneEnforcer(), 'left':MonotoneEnforcer(),'right':MonotoneEnforcer()}
    tracker_index =  [np.where(looper.tr.frames == frame)[0][0] for frame in frame_list]
    times=nar(looper.tr.times[ tracker_index] )
    all_times=looper.tr.times

    #get all the miniscrubbers at once.
    #We should speed this code up.
    #

    if mini_scrubbers is None:
        mini_scrubbers = {}
        for core_id in core_list:
            if mean_velocity:
                do_velocity=True
            else:
                do_velocity=False
            mini_scrubbers[core_id]=  trackage.mini_scrubber(looper.tr,core_id, do_velocity=do_velocity)


    #
    #Loop over all cores and get the bounding box.
    #

    camera.run(core_list, frame_list, mini_scrubbers)



    for frame in frame_list:
        import find_other_cores
        reload(find_other_cores)
        frame_ind = np.where( looper.tr.frames == frame)[0][0]
        left = camera.all_left[frame]
        right = camera.all_right[frame]
        center=camera.all_center[frame]
        position_dict=camera.all_positions[frame]
        big_list, shifts = find_other_cores.centroid_in_box( looper, left, right, mini_scrubbers=mini_scrubbers, last_only=True)
        color_dict = colors.make_core_cmap(big_list, cmap = 'winter', seed = -1)
        color_dict[main_core]=[1.0,0.0,0.0]



        # Check to see if the image was made already,
        # and skips it if it has.
        #outname = "%s/%s_c%04d_n%04d_"%(looper.plot_directory,looper.out_prefix, core_list[0],frame)
        outname = "%s/%s_c%04d_n%04d_xyz"%(looper.plot_directory,looper.out_prefix, main_core,frame)

        ds = looper.load(frame)
        left = camera.all_left[frame]
        right = camera.all_right[frame]
        center=camera.all_center[frame]
        position_dict=camera.all_positions[frame]

        nzones = np.round( (right-left)*2048)
        cg = ds.covering_grid(4, left, nzones)
        reg = ds.region(center,left,right)
        plt.close('all')
        fig,axes=plt.subplots(1,len(axis_list), figsize=(12,4))
        if len(axis_list)==1:
            axes=[axes]
            fig.set_size_inches(4,4)
        for nax,axdir in enumerate(axis_list):
            ax=axes[nax]
            #proj=yt.ProjectionPlot( ds, axdir, field, data_source=reg, center=center, origin='window', weight_field=weight_field)
            #proj.zoom( 1/(right-left).max())
            #proj.save('plots_to_sort/compare')
            Hcoord=ds.coordinates.x_axis[axdir]
            Vcoord=ds.coordinates.y_axis[axdir]
            FFF = cg[field].sum(axdir)
            X = cg['gas','xyz'[Hcoord]]
            Y = cg['gas','xyz'[Vcoord]]
            if axdir in [0,2]:
                FFF = FFF.transpose()
            xlim=[X.min(), X.max()]
            ylim=[Y.min(), Y.max()]
            norm = mpl.colors.LogNorm( FFF.min(), FFF.max())
            ax.imshow( FFF, cmap='Greys',norm=norm, origin='lower',interpolation='nearest',extent=xlim+ylim)
            ax.set(xlabel='xyz'[Hcoord],ylabel='xyz'[Vcoord])

            text_height=0.04
            start_y = 1-text_height
            start_x = 1-4*text_height
            for nc,core_id in enumerate(big_list):
                text_y = start_y - nc*text_height
                text_x = start_x
                points_axis=(text_x,text_y)
                axis_to_data = ax.transAxes + ax.transData.inverted()
                points_data = axis_to_data.transform(points_axis)
                ax.text( text_x,text_y,r'$c%04d$'%core_id, transform=ax.transAxes, color=color_dict.get(core_id,'k'))
                ms = mini_scrubbers[core_id]
                P = np.stack([ms.mean_x[frame_ind], ms.mean_y[frame_ind], ms.mean_z[frame_ind]])
                P += shifts[nc]

                #ax.scatter( P[Hcoord], P[Vcoord])
                ax.plot( [P[Hcoord], points_data[0]], [P[Vcoord], points_data[1]],linewidth=0.1, color=color_dict.get(core_id,'k'))








        fig.tight_layout()
        fig.savefig('plots_to_sort/annotated_%s_c%04d_%s.pdf'%(looper.sim_name,main_core,'xyz'))




