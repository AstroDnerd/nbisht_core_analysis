#
# run get_mountain_tops first
#
from starter2 import *
import xtra_energy
import data_locations as dl
reload(dl)
reload(looper)

if 'this_simname' not in dir():
    this_simname = 'u301'

mountain_top_fname = "datasets_small/%s_mountain_tops_take_9.h5"%this_simname
outname = 'u301_new_tracks_take_9b.h5'


#TO CHANGE
#yt derived field gets defined here
def add_energies(obj):
    if obj.parameters['SelfGravity']:
        def grav_energy(field,data):
            gx =data.ds.arr(grad(data,'PotentialField',0),'code_length/code_time**2')
            gy =data.ds.arr(grad(data,'PotentialField',1),'code_length/code_time**2')
            gz =data.ds.arr(grad(data,'PotentialField',2),'code_length/code_time**2')
            alpha = 1./data.ds['GravitationalConstant'] #=1/4 pi G
            alpha = data.ds.quan(alpha, '1/(code_length**3/code_time**2/code_mass)')
            return ( -(gx**2+gy**2+gz**2)*alpha )
        obj.add_field('grav_energy',grav_energy,validators=[yt.ValidateGridType()],
                 units='code_mass*code_length**2/(code_time**2*code_length**3)', sampling_type='cell')
        def abs_grav_energy(field,data):
            return np.abs( data['grav_energy'] )
        obj.add_field('abs_grav_energy',abs_grav_energy,validators=[yt.ValidateSpatial(1,'PotentialField')],
                 units='code_mass*code_length**2/(code_time**2*code_length**3)', sampling_type='cell')

    def therm_energy(field,data):
        sound_speed = data.ds.quan(1.,'code_velocity')
        rho_0 = data.ds.quan(1.,'code_mass/code_length**2')
        e = sound_speed**2*np.log( data['density']/rho_0)
        therme = data['density']*e
        return therme
    obj.add_field('therm_energy',therm_energy,
                 units='code_mass*code_length**2/(code_time**2*code_length**3)', sampling_type='cell')

if 1:
    """this set of parameters extracts all primitive quantities"""
    target_frame = dl.target_frames[this_simname]
    frame_list =list(range(0,target_frame,10))+[target_frame]
    fields = ['x','y','z','density','velocity_magnitude','magnetic_field_strength', 'velocity_divergence']
    fields += ['velocity_x','velocity_y','velocity_z']
    fields += ['magnetic_field_%s'%s for s in 'xyz']
    fields += ['PotentialField']
    fields += ['particle_pos_x', 'particle_pos_y', 'particle_pos_z', 'particle_index']
    derived=[add_energies] #TO CHANGE

disk_or_new = 'new'
if ('new_looper' not in dir() and disk_or_new) or True:
    new_looper = looper.core_looper(directory= dl.sims[this_simname],
                                     sim_name = this_simname,
                                     out_prefix = this_simname,
                                     target_frame = dl.target_frames[this_simname],
                                     frame_list = frame_list,
                                     core_list =  None,
                                     fields_from_grid=fields,
                                     derived = derived,
                                    do_shift=False
                                  )
    new_looper.plot_directory = "./plots_to_sort"
    new_looper.read_targets(mountain_top_fname)


if 0:
    #verify and remove bad particles.
    #Don't need to do this every time, it takes a minute.
    #need to have the bad_particles file 
    for frame in range(target_frame):
        new_looper.verify_all_particles(frame)
    new_looper.save_bad_particles('%s_bad_particles_full.h5'%this_simname)

if 1:
    new_looper.read_bad_particles('datasets_small/%s_bad_particles_full.h5'%this_simname)
    bad_particle_id = [ 724134,  635702,  661226,  743270,  751995,  718196, 1354060,
                                1362500,  610123,  610189, 1930558, 1046537, 1841352, 1844125,
                                1845574, 1849410, 1853445, 1300291]
    bad_core_id = [ 24,  85,  85,  85,  85,  87, 180, 180, 253, 254, 265, 269, 269,
                                269, 269, 269, 269, 277]
    for bad_core, bad_part in zip(bad_core_id, bad_particle_id):
        new_looper.bad_particles[bad_core]=np.append(
                new_looper.bad_particles[bad_core], bad_part)
    new_looper.remove_bad_particles()
    new_looper.get_tracks()

#this hasn't been tested
def cut_cores_tenfour(loop):
    mountain_top_fname = "datasets_small/%s_mountain_tops_take_9.h5"%mountain_sim
    loops[new_sim].read_targets(mountain_top_fname)

    peaks = nar([loop.targets[core_id].peak_density for core_id in loop.core_list])
    remove = peaks < 1e4
    keep = peaks >= 1e4
    cores_to_cut = loop.core_list[ remove]

    for core_id in cores_to_cut:
        mask = np.where( loop.tr.core_ids == core_id)
        loop.tr.core_ids = np.delete(loop.tr.core_ids, mask)
        loop.tr.particle_ids = np.delete(loop.tr.particle_ids,mask)
        for field in loop.tr.track_dict.keys():
            arr = loop.tr.track_dict[field]
            smaller = np.delete( arr, mask,axis=0)
            loop.tr.track_dict[field]=smaller
    loop.core_list = loop.core_list[ keep]
    for core_id in cores_to_cut:
        del loop.target_indices[core_id]


cut_cores_tenfour( new_looper )


if 1:
    import tracks_read_write
    tracks_read_write.save_loop_trackage_only( new_looper, outname)


if 0:
    print("ZEROS: %d"%( (new_looper.tr.track_dict['density']==0).sum() ))
if 'new_looper' not in dir() and disk_or_new == 'disk':
    print("reload")
    file_list=['TEMP.h5']
    new_looper=looper.core_looper(directory=dl.sims[this_simname])
    new_looper.plot_directory = "/home/dccollins/PigPen"
    for nfile,fname in enumerate(file_list):
        new_looper.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    new_looper.tr.sort_time()
    #import loop_tools
    #loop_tools.re_shift_snaps(new_looper)


if 0:
    import colors
    reload(colors)
    core_cmap = colors.make_core_cmap( new_looper.core_list)
    reload(loop_apps)
    loop_apps.core_proj_multiple(new_looper,#core_list=[258],
                                 #frame_list=[0],
                                 axis_list=[1], color_dict=core_cmap, particles=True,
                                 only_sphere=False,zoom=False,
                                 center_on_sphere=False, annotate=True,tracker_positions=True, shifted_tracker=False)
