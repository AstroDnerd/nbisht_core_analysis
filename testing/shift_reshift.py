from starter2 import *
import xtra_energy
import data_locations as dl
import loop_tools
reload(dl)

reload(looper)
reload(loop_apps)
reload(loop_tools)

this_simname = 'u203'

output_base = "%s_cores"%this_simname
#Cores that overlap more than 90% with u203 core 84
core_list =  np.array([112, 173, 130, 210, 113, 208, 147,  90, 174,  89, 121, 178, 124, 161, 207, 206, 156, 163, 182, 162, 157, 110, 166, 105, 164, 158, 109, 165, 106, 160, 107, 108, 159])
#core_list = [84, 112]
frame_list = [0,10,dl.target_frames[this_simname]] #list(range(0,100,10))+list(range(100,dl.target_frames[this_simname]+1))
frame_list = list(range(0,100,10))+list(range(100,dl.target_frames[this_simname]+1))
#frame_list = range(108)
fields = ['x','y','z','density', 'cell_volume']
derived = []
NewOrDiskOrNeither = 'new'
#tempname = 'u203_core_84_108frames.h5'
tempname = 'u203_core_84_tmp_full.h5'
if 'this_looper' not in dir() and NewOrDiskOrNeither == 'new':
    #for making
    this_looper = looper.core_looper(directory= dl.sims[this_simname],
                                     sim_name = this_simname,
                                     out_prefix = this_simname,
                                     target_frame = dl.target_frames[this_simname],
                                     frame_list = frame_list,
                                     core_list =  core_list,
                                     fields_from_grid=fields,
                                     derived = derived,
                                     do_shift=False
                                  )
    this_looper.plot_directory = "/home/dccollins/PigPen"
    this_looper.get_target_indices(h5_name=dl.peak_list[this_simname])
    #                                 bad_particle_list=dl.bad_particles.get(this_simname,None))
    this_looper.get_tracks()
    #this_looper.make_snapshots()
    if 1:
        this_looper.save(tempname)

if 'this_looper' not in dir() and NewOrDiskOrNeither == 'disk':
    print("reload")
    file_list=[tempname]
    this_looper=looper.core_looper(directory=dl.sims['u05'])
    this_looper.plot_directory = "/home/dccollins/PigPen"
    for nfile,fname in enumerate(file_list):
        this_looper.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    thtr = this_looper.tr
    this_looper.out_prefix='core_84'
    thtr.sort_time()
this_looper.plot_directory = "/home/dccollins/PigPen"

reload(trackage)
def shift_snaps(loop):
    #
    #  Periodic shift of particles is done relative to its final loation.
    #  This is done in trackage.mini_scrubber.scrub.
    #  Here we use the reuslt from the mini scrubber to update the particle positions.
    #
    #  ms: the miniscrubber that takes care of center-of-mass shifts 
    #  ms.this_x, ms.this_y, ms.this_z: shifted GRID positions
    #  snap.pos:  Initially unshifted PARTICLE positions.  These get updated.
    #  pos2:  snap.pos arranged in a single array for easy operations.
    #

    for frame in loop.snaps:
        for core_id in loop.snaps[frame]:

            #Get the objects
            snap = loop.snaps[frame][core_id]
            snap.ds = loop.load(frame)
            tr = loop.tr
            ms = trackage.mini_scrubber( tr, core_id, do_velocity=False)
            ms_index = np.where( tr.frames == frame )[0]
            particle_ids = loop.tr.c([core_id],'particle_id')
            if hasattr(particle_ids,'v'):
                particle_ids = particle_ids.v
            total_particle_index_error = np.abs( particle_ids  - snap.ind.v).sum()
            if total_particle_index_error > 0:
                raise("SORT ERROR")
            pos2 = np.concatenate([ms.this_x[:,ms_index], ms.this_y[:,ms_index], ms.this_z[:,ms_index]], axis=1)

            #The difference between the particle positions (snap.pos) and the shifted grid positions (pos2) 
            #gives the direction of the shift.
            maxshift = np.abs( snap.pos.v - pos2).max() 
            if maxshift > 0.1:
                shift =  pos2 - snap.pos.v
                shift_amount = snap.ds.arr(np.sign(shift), 'code_length')
                to_shift = np.abs( shift) > 0.25
                #if to_shift.sum():
                #    pdb.set_trace()
                snap.pos[ to_shift ] += shift_amount[ to_shift ]

            #Compute the new centroid.
            density = snap.field_values['density']
            volume = snap.field_values['cell_volume']
            m = (density*volume).sum()
            centroid_tmp =  np.array([(snap.pos[:,dim]*density*volume).sum()/m for dim in range(3)])
            snap.R_centroid = snap.ds.arr(centroid_tmp,'cm')


print("Reshift")
shift_snaps(this_looper)
print("done")

if 'color_dict' not in dir():
    color_dict={}
    for core_id in core_list:
        color_dict[core_id] =  np.random.random(3)
#color_dict={84:'r',112:'g'}
#loop_apps.core_proj_follow(this_looper,axis_list=[0], field='PotentialField',
#                           zoom=False, grids=False, particles=False, moving_center=True) #, frame_list=[1])
#slab = {'zmin':zmin, 'zmax':zmax}
#loop_apps.core_proj_follow(this_looper,axis_list=[2], field='PotentialField', slab=slab, zoom=False, only_sphere=False, center_on_sphere=False)#, frame_list=[31])
##                           #zoom=False, grids=False, particles=False, moving_center=True) #, frame_list=[1])
loop_apps.core_proj_multiple(this_looper,axis_list=[2], field='density', color_dict=color_dict, particles=True, fields=False )#, frame_list=[0], moving_center=False)
