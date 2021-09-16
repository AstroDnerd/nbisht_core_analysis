"""
Applications of looper tools.

Decorators:  The @looper.frame_loop and @looper.core_loop decorate functions.
@looper.frame_loop
def this_code(this_looper,other_args):
    ...
will cause this_code to be executed on each frame in this_looper.frame_list.

@looper.core_loop   
def this_code(this_looper):
    stuff

will cause the this_code to be executed for each snapshot on each frame.
"""
from starter2 import *
reload(looper)


def cbump(value):
    """for bumping the centroid value to between 0,1.
    """
    output = (value+0).v
    output[output>1] = output[output>1]-1
    output[output<0] = output[output<0]+1
    return output
@looper.core_loop   #the decorator takes care of the loop
def print_centroid(looper,snapshot): #here's a function, needs to take at least these arguments)
    print("Core %d frame %d centroid (%s)"%(
          snapshot.core_id, snapshot.frame, str( snapshot.R_centroid)))

@looper.frame_loop
def proj_cores(self, axis_list=[0,1,2],core_list=[], field='density'):
    """Full projections of the data, with core particles marked."""
    for axis in axis_list:
        center = self.ds.arr(nar([0.5]*3),'code_length')
        self.proj = self.ds.proj(field,axis,center=center)
        self.proj = yt.ProjectionPlot(self.ds, axis=axis, fields=[field], center=center)
        radius_from_core = []
        core_label = ""
        for nc,core_number in enumerate(core_list):
            core_label += "c%04d_"%core_number
            self.proj.annotate_select_particles(1.0, col='r', indices=self.target_indices[core_number])
            outname = '%s_full_particles_%sn%04d'%(self.out_prefix,core_label,self.current_frame)
        print( self.proj.save(outname))

import annotate_particles_3
reload(annotate_particles_3)
<<<<<<< HEAD
def core_proj_multiple(looper, field='density', axis_list=[0,1,2], color_dict={},
                       force_log=None,linthresh=100,
                       core_list=None,frame_list=None, clobber=True,
                       grids=True, particles=True, 
                       fields=False, velocity=False, code_length=True, lic=False, annotate=False, 
                       moving_center=False, only_sphere=True, center_on_sphere=True, slab=None, zoom=True, 
                       marker_size = 1,
                      tracker_positions=True, shifted_tracker=True, float_positions=False):
=======
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
    def __init__(self):
        self.values=[]
        self.extrema=[]
        self.sign=None
    def __call__(self,val,verbose=False):
        if len(self.values) < 2:
            new_val = val
        else:
            delta_i = val - self.values[-1]
            delta_im1 = self.values[-1]-self.values[-2]
            if self.sign is None:
                if (delta_im1 ==0).all():
                    if verbose:
                        print("I'm not sure what to do if delta_im1==0.  Check that the colorbar is right.")
                    delta_im1 = delta_i +0
                self.sign = np.sign(delta_im1)
            if len(self.values) > 3:
                arr = np.column_stack(self.values[-4:])
                this_sign = np.zeros(len(val))
                for nv in np.arange( len(val)):
                    pfit = np.polyfit( np.arange(4), arr[nv,:],1)
                    this_sign=np.sign(pfit[0])
                self.sign=this_sign
                if verbose:
                    print(self.sign)
                



            direction = np.sign( delta_i * self.sign)
            new_val = copy.copy(val)
            new_val[ direction <= 0] = self.extrema[-1][direction<=0]
            #if verbose:
            #    print( "MONO pre ", self.prev[-1])
            #    print( "MONO val ", val)
            #    print( "MONO val ", new_val)
            #    print( "MONO pre ", self.prev)
            #    print( "MONO dir ", direction)
            #    print( "MONO delm", delta_im1)
            #    print( "MONO deli", delta_i)
        self.extrema.append(new_val)
        self.values.append(val)
        return new_val




def annoying_limiter(looper, field='density', axis_list=[0,1,2], color_dict={},force_log=None,linthresh=100,
                    core_list=None,frame_list=None, clobber=True,
                       only_sphere=True, center_on_sphere=True,
                       zoom=True, moving_center=False, slab=None,
                       grids=True, particles=True, annotate=False, 
                      fields=False, velocity=False, lic=False, 
                       code_length=True, 
                      tracker_positions=True, shifted_tracker=True, monotonic=False, float_positions=False,
                      marker_size=1):
>>>>>>> 40e1c991ca757f3dced34474fd2958b18fc829ea
    if core_list is None:
        core_list = looper.core_list
    if frame_list is None:
        frame_list = looper.frame_list
    all_png = glob.glob("*png")
    tr = looper.tr
    if monotonic is True:
        monotonic  ={'zlim':MonotoneEnforcer2(), 'left':MonotoneEnforcer2(),'right':MonotoneEnforcer2()}
        #monotonic  ={'zlim':MonotoneEnforcer(), 'left':MonotoneEnforcer(),'right':MonotoneEnforcer()}
    mini_scrubbers = {}
    for core_id in core_list:
        mini_scrubbers[core_id]=  trackage.mini_scrubber(looper.tr,core_id, do_velocity=False)
    for frame in frame_list:
        ds = looper.load(frame)

        # Check to see if the image was made already,
        # and skips it if it has.
        outname = "%s/%s_n%04d_multi"%(looper.plot_directory,looper.out_prefix, frame)
        got_one = False
        for i in all_png:
            if i.startswith(outname):
                print(i)
                got_one=True
        if got_one and not clobber:
            print("File exists, skipping")
            return
        #get extents and bounding region
<<<<<<< HEAD
        left =  np.array([1,1,1])
        right = np.array([0,0,0])
=======
        #it's backwards because we're looking for extrema
        left =  ds.domain_right_edge.v
        right = ds.domain_left_edge.v

        #
        # Find the extents of all cores.
        # Fill position array.
        #
>>>>>>> 40e1c991ca757f3dced34474fd2958b18fc829ea
        position_dict={}
        for core_id in core_list:
            ms = mini_scrubbers[core_id]
            if tracker_positions:
                frame_ind = np.where(looper.tr.frames == frame)[0]
                if shifted_tracker:
                    this_x=ds.arr(ms.this_x[:,frame_ind],"code_length")
                    this_y=ds.arr(ms.this_y[:,frame_ind],"code_length")
                    this_z=ds.arr(ms.this_z[:,frame_ind],"code_length")
                else:
                    this_x=ds.arr(ms.raw_x[:,frame_ind],"code_length")
                    this_y=ds.arr(ms.raw_y[:,frame_ind],"code_length")
                    this_z=ds.arr(ms.raw_z[:,frame_ind],"code_length")
                if float_positions:
                    ms.make_floats(core_id)
                    this_x = ds.arr( ms.float_x[:,frame_ind], 'code_length')
                    this_y = ds.arr( ms.float_y[:,frame_ind], 'code_length')
                    this_z = ds.arr( ms.float_z[:,frame_ind], 'code_length')

                positions = np.column_stack([this_x,this_y,this_z])
                position_dict[core_id] = positions
            else:
                snapshot = looper.snaps[frame][core_id]
                if snapshot.R_centroid is None:
                    snapshot.get_all_properties()
                positions = snapshot.pos
                position_dict[core_id]=positions
            this_left =  positions.min(axis=0)
            this_right = positions.max(axis=0)
            left = np.row_stack([this_left,left]).min(axis=0)
            right = np.row_stack([this_right,right]).max(axis=0)

        center = 0.5*(left+right)
        left = np.row_stack([left,center - 1/128]).min(axis=0)
        right = np.row_stack([right,center + 1/128]).max(axis=0)
        left = ds.arr(left,'code_length')
        right = ds.arr(right,'code_length')
        center = 0.5*(left+right)
        center = ds.arr(center,'code_length')


        if not center_on_sphere:
            center = 0.5*(ds.domain_left_edge+ds.domain_right_edge)

        if monotonic:
            left = monotonic['left'](left, verbose=True)
            right  = monotonic['right'](right)
            center = 0.5*(left+right)

            #jleft = np.row_stack([this_left,left]).min(axis=0)
            #jright = np.row_stack([this_right,right]).max(axis=0)
    fig,ax=plt.subplots(1,1, figsize=(12,12))
    for core_id in core_list:
        for ny in np.arange(mini_scrubbers[core_id].nparticles):
            ax.plot(mini_scrubbers[core_id].float_y[ny,:],linewidth=0.1)
        RRR=np.column_stack(monotonic['right'].values)
        LLL=np.column_stack(monotonic['left'].values)
        ax.plot( RRR[1,:],c='r')
        ax.plot( LLL[1,:],c='k')

    fig.savefig('plots_to_sort/Y.png')
    return monotonic, RRR


        #
        # main plot loop



def core_proj_multiple(looper, field='density', axis_list=[0,1,2], color_dict={},force_log=None,linthresh=100,
                    core_list=None,frame_list=None, clobber=True,
                       only_sphere=True, center_on_sphere=True,
                       zoom=True, moving_center=False, slab=None,
                       grids=True, particles=True, annotate=False, 
                      fields=False, velocity=False, lic=False, 
                       code_length=True, 
                      tracker_positions=True, shifted_tracker=True, monotonic=False, float_positions=False,
                      marker_size=1):
    if core_list is None:
        core_list = looper.core_list
    if frame_list is None:
        frame_list = looper.frame_list
    all_png = glob.glob("*png")
    tr = looper.tr
    if monotonic is True:
        monotonic  ={'zlim':MonotoneEnforcer2(), 'left':MonotoneEnforcer2(),'right':MonotoneEnforcer2()}
        #monotonic  ={'zlim':MonotoneEnforcer(), 'left':MonotoneEnforcer(),'right':MonotoneEnforcer()}
    mini_scrubbers = {}
    fig,ax=plt.subplots(1,1)
    for core_id in core_list:
        mini_scrubbers[core_id]=  trackage.mini_scrubber(looper.tr,core_id, do_velocity=False)
        for ny in np.arange(mini_scrubbers[core_id].nparticles):
            ax.plot(mini_scrubbers[core_id].this_y[ny,:])
    fig.savefig('plots_to_sort/Y.png')
    for frame in frame_list:
        ds = looper.load(frame)

        # Check to see if the image was made already,
        # and skips it if it has.
        outname = "%s/%s_n%04d_multi"%(looper.plot_directory,looper.out_prefix, frame)
        got_one = False
        for i in all_png:
            if i.startswith(outname):
                print(i)
                got_one=True
        if got_one and not clobber:
            print("File exists, skipping")
            return
        #get extents and bounding region
        #it's backwards because we're looking for extrema
        left =  ds.domain_right_edge.v
        right = ds.domain_left_edge.v

        #
        # Find the extents of all cores.
        # Fill position array.
        #
        position_dict={}
        for core_id in core_list:
            ms = mini_scrubbers[core_id]
            if tracker_positions:
                frame_ind = np.where(looper.tr.frames == frame)[0]
                if shifted_tracker:
                    this_x=ds.arr(ms.this_x[:,frame_ind],"code_length")
                    this_y=ds.arr(ms.this_y[:,frame_ind],"code_length")
                    this_z=ds.arr(ms.this_z[:,frame_ind],"code_length")
                else:
                    this_x=ds.arr(ms.raw_x[:,frame_ind],"code_length")
                    this_y=ds.arr(ms.raw_y[:,frame_ind],"code_length")
                    this_z=ds.arr(ms.raw_z[:,frame_ind],"code_length")
                if float_positions:
                    ms.make_floats(core_id)
                    this_x = ds.arr( ms.float_x[:,frame_ind], 'code_length')
                    this_y = ds.arr( ms.float_y[:,frame_ind], 'code_length')
                    this_z = ds.arr( ms.float_z[:,frame_ind], 'code_length')

                positions = np.column_stack([this_x,this_y,this_z])
                position_dict[core_id] = positions
            else:
                snapshot = looper.snaps[frame][core_id]
                if snapshot.R_centroid is None:
                    snapshot.get_all_properties()
                positions = snapshot.pos
            this_left =  positions.min(axis=0)
            this_right = positions.max(axis=0)
            left = np.row_stack([this_left,left]).min(axis=0)
            right = np.row_stack([this_right,right]).max(axis=0)

        center = 0.5*(left+right)
        left = np.row_stack([left,center - 1/128]).min(axis=0)
        right = np.row_stack([right,center + 1/128]).max(axis=0)
        left = ds.arr(left,'code_length')
        right = ds.arr(right,'code_length')
        center = 0.5*(left+right)
        center = ds.arr(center,'code_length')

        if not center_on_sphere:
            center = 0.5*(ds.domain_left_edge+ds.domain_right_edge)

<<<<<<< HEAD
=======
        if monotonic:
            left = monotonic['left'](left)
            right  = monotonic['right'](right)
            center = 0.5*(left+right)

            #jleft = np.row_stack([this_left,left]).min(axis=0)
            #jright = np.row_stack([this_right,right]).max(axis=0)


        #
        # main plot loop
        #
>>>>>>> 40e1c991ca757f3dced34474fd2958b18fc829ea
        for ax in axis_list:
            Rmax = np.sqrt( ( (right-left)**2).sum(axis=0)).max()
            scale = Rmax.v #2*max([Rmax, scale_min]).v
            print("SCALE", scale)
            scale = min([scale,1])
            if only_sphere:
                sph = ds.region(center,left,right)
                #sph = ds.sphere(center,Rmax)
                proj = ds.proj(field,ax,center=center, data_source = sph) 
            else:
                proj = ds.proj(field,ax,center=center)
            looper.proj=proj
            #if moving_center:
            #    #not quite working.
            #    pw = proj.to_pw(center=[0.5]*3, origin = 'native')#center = center,width=(1.0,'code_length'))
            #    pw.set_center([center[0],center[1]])
            #else:
            #    pw = proj.to_pw(center = center,width=(1.0,'code_length'), origin='domain')
            pw = proj.to_pw(center = center, origin='domain')
            
            if zoom:
<<<<<<< HEAD
                pw.zoom(1./(2*scale))

            pw.set_cmap(field,'gray')
=======
                pw.zoom(1./(scale))
            if monotonic:
                array = proj[field]
                used_min = proj[field][ proj['weight_field']>0].min()
                zmin,zmax = used_min, proj[field].max()
                zlim = nar([zmin,zmax])
                zlim = monotonic['zlim'](zlim, verbose=True)
                pw.set_zlim(field,zlim[0],zlim[1])
                print("ZZZZ %0.2e %0.2e MMM %0.2e %0.2e"%(zmin,zmax,zlim[0],zlim[1]))

            pw.set_cmap(field,'Greys')
>>>>>>> 40e1c991ca757f3dced34474fd2958b18fc829ea
            if force_log is not None:
                pw.set_log(field,force_log,linthresh=linthresh)
            for core_id in core_list:
                positions = position_dict[core_id]
                color=color_dict[core_id]
                if annotate:
                    #pw.annotate_sphere(snapshot.R_centroid,Rmax, circle_args={'color':color} ) #R_mag.max())
                    centroid=positions.mean(axis=0).v
                    pw.annotate_text(centroid,
                                     "%d"%core_id,text_args={'color':color}, 
                                     inset_box_args={'visible':False},
                                     coord_system='data')
                if particles:
                    pw.annotate_these_particles2(1.0, col=[color]*positions.shape[0], positions=positions, 
                                                 p_size=marker_size)
        pw.save('plots_to_sort/P3')
        if lic:
            pw.annotate_line_integral_convolution('magnetic_field_x','magnetic_field_y', lim=(0.5,0.65))
        if fields:
            #pw.annotate_magnetic_field(plot_args={'color':'b'})
            pw.annotate_streamlines("magnetic_field_x","magnetic_field_y",plot_args={'color':'b'})
        if velocity:
            pw.annotate_velocity()
        if grids:
            pw.annotate_grids()
        if code_length:
            pw.set_axes_unit('code_length')
        if field in ['MagVelAngle']:
            pw.set_cmap('MagVelAngle','hsv')
            pw.set_zlim('MagVelAngle',0,180)
        looper.pw = pw
        print(pw.save(outname))


#
# Other functions are not as useful.
#
import annotate_particles_3
reload(annotate_particles_3)
from scipy.spatial import ConvexHull
@looper.frame_loop
def proj_cores2(self, axis_list=[0,1,2],core_list=[], field='density',color='r'):
    """Full projections of the data, with core particles marked."""
    for snapshot in self.snaps[self.current_frame].values():
        if snapshot.R_centroid is None:
            snapshot.get_all_properties()
    for axis in axis_list:
        ds = self.ds_list[self.current_frame]
        proj_center = ds.arr([0.12768555, 0.4987793 , 0.17797852], 'code_length')
        proj_center = ds.arr([0.6,0.6,0.6], 'code_length')

        center = self.ds.arr(nar([0.5]*3),'code_length')
        print("poot")
        self.proj = yt.ProjectionPlot(self.ds, axis=axis, fields=[field])#, center=center)
        self.proj.set_cmap(field,'Greys')
        radius_from_core = [] 
        core_label = ""
        for nc,core_number in enumerate(core_list):
            ms = trackage.mini_scrubber(self.tr,core_number)
            frame_ind = np.where(self.tr.frames == self.current_frame)[0]
            this_x=ds.arr(ms.raw_x[:,frame_ind],"code_length")
            this_y=ds.arr(ms.raw_y[:,frame_ind],"code_length")
            this_z=ds.arr(ms.raw_z[:,frame_ind],"code_length")
            points3d = this_x,this_y,this_z

            snapshot = self.snaps[self.current_frame][core_number]
            center = ds.arr(snapshot.R_centroid,'code_length')
            print("Core %d center %s"%(core_number, str(center)))
            core_label += "c%04d_"%core_number
            self.proj.annotate_select_particles4(1.0, col='r', indices=self.target_indices[core_number])
            outname = '%s_full_particles_%sn%04d'%(self.out_prefix,core_label,self.current_frame)
            reload( annotate_particles_3)
            self.proj.annotate_text(center,
                             "%d"%snapshot.core_id,text_args={'color':color}, 
                             inset_box_args={'visible':False},
                             coord_system='data')
            self.proj.annotate_convex_hull_1(1,points=points3d)
        print( self.proj.save(outname))

def proj_select_particles(this_looper, frame_list=None, axis_list=[0,1,2],core_list=[], field='density',color='r'):
    """Full projections of the data, with core particles marked."""
    rm = rainbow_map(len(core_list))
    if frame_list is None:
        frame_list = this_looper.frame_list

    for frame in frame_list:

        for snapshot in this_looper.snaps[frame].values():
            if snapshot.R_centroid is None:
                snapshot.get_all_properties()
        for axis in axis_list:
            ds = this_looper.load(frame)
            proj_center = ds.arr([0.12768555, 0.4987793 , 0.17797852], 'code_length')
            proj_center = ds.arr([0.6,0.6,0.6], 'code_length')

            center = this_looper.ds.arr(nar([0.5]*3),'code_length')
            print("poot")
            this_looper.proj = yt.ProjectionPlot(this_looper.ds, axis=axis, fields=[field])#, center=center)
            this_looper.proj.set_cmap(field,'Greys')
            core_label = ""
            mean_center=0
            for nc,core_number in enumerate(core_list):
                c=rm(nc)
                snapshot = this_looper.snaps[frame][core_number]
                center = ds.arr(snapshot.R_centroid,'code_length')
                mean_center = center + mean_center
                print("Core %d center %s"%(core_number, str(center)))
                core_label += "c%04d_"%core_number
                this_looper.proj.annotate_select_particles4(1.0, indices=this_looper.target_indices[core_number],col=c)
                reload( annotate_particles_3)
                this_looper.proj.annotate_text(center,
                                 "%d"%snapshot.core_id,text_args={'color':c}, 
                                 inset_box_args={'visible':False},
                                 coord_system='data')
            mean_center /= len(core_list)
            plot_x = [1,2,0][axis]
            plot_y = [2,0,1][axis]
            this_center = mean_center[plot_x],mean_center[plot_y]
            this_looper.proj.set_center(this_center)
            this_looper.proj.zoom(4)

            outname = 'plots_to_sort/%s_full_particles_%sn%04d'%(this_looper.out_prefix,"cores_10_11_",frame)
            print( this_looper.proj.save(outname))

@looper.frame_loop
def proj_cores(self, axis_list=[0,1,2],core_list=[], field='density'):
    """Full projections of the data, with core particles marked."""
    for axis in axis_list:
        center = self.ds.arr(nar([0.5]*3),'code_length')
        self.proj = self.ds.proj(field,axis,center=center)
        self.proj = yt.ProjectionPlot(self.ds, axis=axis, fields=[field], center=center)
        radius_from_core = []
        core_label = ""
        for nc,core_number in enumerate(core_list):
            core_label += "c%04d_"%core_number
            self.proj.annotate_select_particles(1.0, col='r', indices=self.target_indices[core_number])
            outname = '%s_full_particles_%sn%04d'%(self.out_prefix,core_label,self.current_frame)
        print( self.proj.save(outname))

@looper.frame_loop
def select_particles(looper,these_particles=None,axis_list=[0,1,2]):
    for axis in axis_list:
        proj = looper.ds.proj('density',axis,center='c')
        pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap('density','gray')
        radius_from_core = []
        core_label = ""
        pw.annotate_select_particles(1.0, col='r', indices=these_particles, p_size=5)
        outname = '%s_select_n%04d'%(looper.out_prefix,looper.current_frame)
        print(outname)
        print( pw.save(outname))

def core_proj_follow_b(looper, field='density', axis_list=[0,1,2], color='r',force_log=None,linthresh=100,
                    core_list=None,frame_list=None, clobber=True,zoom=True, grids=True, particles=True, moving_center=False, 
                    only_sphere=True, center_on_sphere=True, slab=None, annotate=True, p_size=0.1):
    print("THIS CODE HAS BEEN REMOVED IN FAVOR OF core_proj_multiple")


@looper.core_loop
def core_circle(looper,snapshot, axis_list=[0,1,2]):
    for axis in axis_list:
        proj = looper.ds.proj('density',axis,center='c')
        pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap('density','gray')
        radius_from_core = []
        core_label = ""
        pw.annotate_select_particles(1.0, col='r', indices=looper.target_indices[looper.core_id])
        pw.annotate_sphere(snapshot.R_centroid,snapshot.R_mag.max() ) #R_mag.max())
        outname = '%s_fullcircle_cores_%sn%04d'%(looper.out_prefix,looper.core_id,looper.current_frame)
        print( pw.save(outname))

@looper.frame_loop
def proj_cores_with_annotations(self, axis_list=[0,1,2], color_dict={}):
    """Full projections, with particles plotted.  Also plots the radius and number for each core.
    Colors can be passed in through a dictionary, indexed by core_id"""
    for axis in axis_list:
        proj = self.ds.proj('density',axis,center='c')
        pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap('density','gray')
        radius_from_core = []
        core_label = ""
        for nc,core_number in enumerate(self.target_indices):
            this_snapshot = self.make_snapshot(self.current_frame,core_number)
            if this_snapshot.R_centroid is None:
                this_snapshot.get_all_properties()
            center = this_snapshot.R_centroid
            color = color_dict.get(core_number,'r')
            core_label += "c%04d_"%core_number
            pw.annotate_text(cbump(center),
                             "%d"%core_number,text_args={'color':color}, 
                             inset_box_args={'visible':False},
                             coord_system='data')
            pw.annotate_select_particles(1.0, col=color, indices=self.target_indices[core_number])
        outname = '%s_full_particles_%sn%04d'%(self.out_prefix,core_label,self.current_frame)
        print( pw.save(outname))

@looper.frame_loop
def proj_with_species(self, field='density',axis_list=[0,1,2], color_dict={}, core_list=None):
    """Full projections, with particles plotted.  Also plots the radius and number for each core.
    Colors can be passed in through a dictionary, indexed by core_id"""
    if core_list is None:
        core_list = self.target_indices.keys()
    for axis in axis_list:
        proj = self.ds.proj(field,axis,center='c')
        pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap('density','gray')
        for nc,core_number in enumerate(core_list):
            this_snapshot = self.make_snapshot(self.current_frame,core_number)
            if this_snapshot.R_centroid is None:
                this_snapshot.get_all_properties()
            center = this_snapshot.R_centroid
            color = color_dict.get(core_number,'r')
            pw.annotate_text(cbump(center),
                             "%d"%core_number,text_args={'color':color}, 
                             inset_box_args={'visible':False},
                             coord_system='data')
        outname = '%s_proj_regime_n%04d'%(self.out_prefix,self.current_frame)
        print( pw.save(outname))

def phase_with_preimage(this_looper,frame,fields,weight_field=None, bins=None, xlim=None,ylim=None,zlim=None):
    print('yay')
    ds = this_looper.load(frame)
    ad = ds.all_data()
    phase_all = yt.create_profile(ad,fields[:2], fields[2],weight_field=None, override_bins=bins)
    p = yt.ParticlePhasePlot.from_profile(phase_all)
    if xlim:
        p.set_xlim(xlim[0],xlim[1])
    if ylim:
        p.set_ylim(ylim[0],ylim[1])
    if zlim:
        p.set_zlim(fields[2],zlim[0],zlim[1])

    #this save is important.
    p.save("plots_to_sort/%s"%this_looper.out_prefix)
    tr=this_looper.tr
    location = np.where(tr.frames==frame)[0][0]
    the_x = tr.track_dict[fields[0]][:,location]
    the_y = tr.track_dict[fields[1]][:,location]
    cell_plot = p.plots['cell_volume']
    this_axes = cell_plot.axes
    this_axes.scatter(the_x,the_y,c='k',s=0.5)
    print(p.save("plots_to_sort/%s_n%04d"%(this_looper.out_prefix,frame)))




