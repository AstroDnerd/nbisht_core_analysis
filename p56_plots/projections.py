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

#Look for app_test.h5 and read from that.
#If not, make a new looper and write app_test.h5.
if 'this_looper' not in dir():

    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    this_looper = looper.core_looper(directory= directory, 
                                     derived=[xtra_energy.add_force_terms], 
                                     sim_name = 'u05', 
                                     out_prefix = 'plots_to_sort/u05', 
                                     target_frame = 125, 
                                     frame_list = frame_list, 
                                     core_list = core_list, 
                                     fields_from_grid=['x','y','z']+fields)
    this_looper.get_target_indices(h5_name='datasets_small/u05_0125_peaklist.h5',
                                     bad_particle_list='bad_particles.h5')
    this_looper.get_tracks()

@looper.frame_loop
def proj_cores_annotate_zoom(self, axis_list=[0,1,2],core_list=None, center_cores=[],field='density',color='r',
                            zoom_level=4,cb_label=None):
    """Full projections of the data, with core particles marked."""
    if core_list is None:
        core_list = list(self.target_indices.keys())
    for snapshot in self.snaps[self.current_frame].values():
        if snapshot.R_centroid is None:
            snapshot.get_all_properties()
    for axis in axis_list:
        ds = self.ds_list[self.current_frame]
        proj_center = ds.arr([0.6,0.6,0.6], 'code_length')

        center = self.ds.arr(nar([0.5]*3),'code_length')
        print("poot")
        #self.proj = yt.ProjectionPlot(self.ds, axis=axis, fields=[field])#, center=center)
        proj = ds.proj(field,axis,center=center)
        self.proj=proj.to_pw(width=(1.0,'code_length'),center=center,origin='domain')
        #self.proj.set_cmap(field,'Greys')
        self.proj.set_cmap(field,'Greys')
        self.proj.zoom(zoom_level)
        self.proj.set_axes_unit('code_length')
        for nc,core_number in enumerate(core_list):
            ms = trackage.mini_scrubber(self.tr,core_number)
            snapshot = self.snaps[self.current_frame][core_number]
            center = ds.arr(snapshot.R_centroid,'code_length')
            print("Core %d center %s"%(core_number, str(center)))
            #self.proj.annotate_select_particles4(1.0, col='r', indices=self.target_indices[core_number])
            self.proj.annotate_text(loop_apps.cbump(center),
                             "%d"%snapshot.core_id,text_args={'color':color}, 
                             inset_box_args={'visible':False},
                             coord_system='data')

        for core_number in center_cores:
            snapshot = self.snaps[self.current_frame][core_number]
            center = ds.arr(snapshot.R_centroid,'code_length')
            plot_x = [1,2,0][axis]
            plot_y = [2,0,1][axis]
            this_center = center[plot_x],center[plot_y]
            outname = '%s_core_zoom_annotate_c%04d_n%04d'%(self.out_prefix,core_number,self.current_frame)
            if cb_label is not None:
                self.proj.set_colorbar_label(field,cb_label)
            self.proj.set_center(this_center)
            print( self.proj.save(outname))

if 1:

    proj_cores_annotate_zoom(this_looper,axis_list=[0],core_list=core_list,center_cores=[8],cb_label='density')

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

