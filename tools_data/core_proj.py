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

class MonotoneEnforcer():
    def __init__(self):
        self.prev=[]
        self.sign=None
    def __call__(self,val,verbose=False):
        if len(self.prev) < 2:
            new_val = val
        else:
            delta_i = val - self.prev[-1]
            delta_im1 = self.prev[-1]-self.prev[-2]
            if self.sign is None:
                if (delta_im1 ==0).all():
                    if verbose:
                        print("I'm not sure what to do if delta_im1==0.  Check that the colorbar is right.")
                    delta_im1 = delta_i +0
                self.sign = np.sign(delta_im1)
            direction = np.sign( delta_i * self.sign)
            print(self.sign)
            new_val = copy.copy(val)
            new_val[ direction <= 0] = self.prev[-1][direction<=0]
            if verbose:
                print( "MONO pre ", self.prev[-1])
                print( "MONO val ", val)
                print( "MONO val ", new_val)
                print( "MONO pre ", self.prev)
                print( "MONO dir ", direction)
                print( "MONO delm", delta_im1)
                print( "MONO deli", delta_i)
        self.prev.append(new_val)
        return new_val

class MonotoneEnforcer2():
    def __init__(self, nvalues=4):
        self.values=[]
        self.extrema=[]
        self.sign=None
        self.nvalues=nvalues
    def __call__(self,val,verbose=False):
        if len(self.values) < self.nvalues:
            new_val = val
        else:
            if self.sign is None:
                self.sign = np.zeros(len(val))
            delta_i = val - self.extrema[-1]

            if len(self.values) >= self.nvalues:
                arr = np.column_stack(self.values[-self.nvalues:])
                this_sign = np.zeros(len(val))
                for nv in np.arange( len(val)):
                    pfit = np.polyfit( np.arange(self.nvalues), arr[nv,:],1)
                    this_sign=np.sign(pfit[0])
                    self.sign[nv]=this_sign
                if verbose:
                    print(self.sign)
                
            direction = np.sign( delta_i * self.sign)
            new_val = copy.copy(val)
            new_val[ direction <= 0] = self.extrema[-1][direction<=0]
        self.extrema.append(new_val)
        self.values.append(val)
        return new_val

def core_proj_multiple(looper, camera=None, field='density', axis_list=[0,1,2], color_dict={},force_log=None,linthresh=100,
                       main_core=None,
                    core_list=None,frame_list=None, clobber=True, weight_field=None,
                       only_sphere=True, center_on_sphere=True,
                       zoom=True, moving_center=False, slab=None,
                       grids=True, plot_particles=True, annotate=False, 
                      fields=False, velocity=False, mean_velocity=False, lic=False, 
                       code_length=True, 
                      tracker_positions=True, shifted_tracker=True, monotonic=False, float_positions=False,
                      marker_size=1, path_only=False, verbose=False, plot_y_tracks=True, plot_points=False,
                       zlim=None,
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

    multiproj = False
    if len(axis_list):
        multiproj=True


    if not clobber:
        work_to_do=False
        things=[]
        
        for frame in frame_list:

            # Check to see if the image was made already,
            # and skips it if it has.
            #outname = "%s/%s_c%04d_n%04d_"%(looper.plot_directory,looper.out_prefix, core_list[0],frame)
            outname = "%s/%s_c%04d_n%04d_xyz"%(looper.plot_directory,looper.out_prefix, main_core,frame)
            got_this=False
            for i in all_png:
                if i.startswith(outname):
                    got_this=True
            if not got_this:
                work_to_do=True
                things.append(outname)

    if not work_to_do:
        print("Nothing to do, exit")
        return -1
    #
    #get all the miniscrubbers at once.
    #We should speed this code up.
    #

    if verbose:
        print("Mini scrubbers")
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

        # Check to see if the image was made already,
        # and skips it if it has.
        #outname = "%s/%s_c%04d_n%04d_"%(looper.plot_directory,looper.out_prefix, core_list[0],frame)
        outname = "%s/%s_c%04d_n%04d_xyz"%(looper.plot_directory,looper.out_prefix, main_core,frame)
        got_one = False
        for i in all_png:
            if i.startswith(outname):
                print(i)
                got_one=True
        if got_one and not clobber:
            print("File exists, skipping")
            continue
        ds = looper.load(frame)
        left = camera.all_left[frame]
        right = camera.all_right[frame]
        center=camera.all_center[frame]
        position_dict=camera.all_positions[frame]

        #
        # main plot loop
        #

        plot_collector=[]
        sph = ds.region(center,left,right)
        for ax in axis_list:
            Rmax = np.sqrt( ( (right-left)**2).max(axis=0)).max()
            if hasattr(Rmax,'v'):
                scale = Rmax.v #2*max([Rmax, scale_min]).v
            else:
                scale = Rmax

            print("SCALE", scale)
            scale = min([scale,1])
            scale = max([scale,4/128])
            if 0:
                if only_sphere:
                    sph = ds.region(center,left,right)
                    #sph = ds.sphere(center,Rmax)
                    proj = ds.proj(field,ax,center=center, data_source = sph, weight_field=weight_field)
                else:
                    bv = ds.arr([10,10,10],'code_velocity')
                    proj = ds.proj(field,ax,center=center, field_parameters={'bulk_velocity':bv})
            if 1:
                if only_sphere:
                    #sph = ds.sphere(center,Rmax)
                    #proj = ds.proj(field,ax,center=center, data_source = sph, weight_field=weight_field)
                    pw = yt.ProjectionPlot(ds, ax, field, data_source=sph, center=center, origin='window', weight_field=weight_field)
                    plot_collector.append(pw)
                else:
                    bv = ds.arr([10,10,10],'code_velocity')
                    proj = ds.proj(field,ax,center=center, field_parameters={'bulk_velocity':bv})
            if zoom == 'scale' or zoom is True:
                pw.zoom(1./(scale))
            elif type(zoom) == float or type(zoom) == int:
                pw.zoom(zoom)
            if monotonic:
                array = proj[field]
                used_min = proj[field][ proj['weight_field']>0].min()
                zmin,zmax = used_min, proj[field].max()
                zlim = nar([zmin,zmax])
                zlim = monotonic['zlim'](zlim, verbose=True)
                pw.set_zlim(field,zlim[0],zlim[1])
                print("ZZZZ %0.2e %0.2e MMM %0.2e %0.2e"%(zmin,zmax,zlim[0],zlim[1]))
            if zlim is not None:
                pw.set_zlim(field,zlim[0],zlim[1])

            pw.set_cmap(field,cmap)
            if force_log is not None:
                pw.set_log(field,force_log,linthresh=linthresh)
            for core_id in core_list:
                positions = position_dict[core_id]
                color=color_dict[core_id]
                #color = [1.0,0.0,0.0,0.8]
                #alpha=0.1
                alpha=0.4
                marker_size=1
                if annotate:
                    #pw.annotate_sphere(snapshot.R_centroid,Rmax, circle_args={'color':color} ) #R_mag.max())
                    centroid=positions.mean(axis=0)
                    pw.annotate_text(centroid,
                                     "%d"%core_id,text_args={'color':color}, 
                                     inset_box_args={'visible':False},
                                     coord_system='data')
                if plot_particles:
                    pw.annotate_these_particles4(1.0, col=[color]*positions.shape[0], positions=positions, 
                                                 p_size=marker_size, alpha=alpha)
            if lic:
                pw.annotate_line_integral_convolution('magnetic_field_x','magnetic_field_y', lim=(0.5,0.65))
            if fields:
                #pw.annotate_magnetic_field(plot_args={'color':'b'})
                pw.annotate_streamlines("magnetic_field_x","magnetic_field_y",plot_args={'color':'b'})
            if velocity:
                pw.annotate_velocity()
            if mean_velocity:
                Hcoord=ds.coordinates.x_axis[ax]
                Vcoord=ds.coordinates.y_axis[ax]
                VH = "velocity_"+"xyz"[Hcoord]
                VV = "velocity_"+"xyz"[Vcoord]
                ms = mini_scrubbers[core_id]
                ms.get_central_at_once(core_id)
                bulk = np.stack([ms.mean_vx,ms.mean_vy,ms.mean_vz])[:,frame_ind]
                #bulk = np.stack([ms.vxc,ms.vyc,ms.vzc])[:,frame_ind]
                def vvh(field,data):
                    thing1=data[VH]
                    thing2=ds.arr(bulk[Hcoord],'code_velocity')
                    return thing1-thing2
                def vvv(field,data):
                    return data[VV]-ds.arr(bulk[Vcoord],'code_velocity')
                ds.add_field(('gas','vh'),vvh,units='code_velocity',sampling_type='cell')
                ds.add_field(('gas','vv'),vvv,units='code_velocity',sampling_type='cell')

                pw.annotate_quiver('vh','vv')


            if grids:
                pw.annotate_grids()
            if code_length:
                pw.set_axes_unit('code_length')
            if field in ['MagVelAngle']:
                pw.set_cmap('MagVelAngle','hsv')
                pw.set_zlim('MagVelAngle',0,180)
            #looper.pw = pw
            if not multiproj:
                raise
        if multiproj:
            from mpl_toolkits.axes_grid1 import AxesGrid
            if multiproj:
                fig=plt.figure()
                grid = AxesGrid(
                fig,(0.00,0.00,1.0,1.0),
                nrows_ncols=(1, 3),
                axes_pad=0.0,
                label_mode="1",
                share_all=True)

                for axdir in axis_list:
                    #plot = p.plots[field]
                    plot = plot_collector[axdir].plots[YT_density]
                    plot.figure = fig
                    plot.axes = grid[axdir].axes
                    plot.cax = grid.cbar_axes[axdir]
                for axdir in axis_list:
                    plot_collector[axdir]._setup_plots()
                for axdir in [0,1,2]:
                    grid[axdir].set( xticks=[],yticks=[], xlabel='',ylabel='')
                fig.set_size_inches(12,4)
            outname = "%s/%s_c%04d_n%04d_xyz"%(looper.plot_directory,looper.out_prefix, main_core,frame)

            print(outname)
            fig.savefig(outname)
            #print(pw.save(outname))
        del sph
        plt.close(fig)

