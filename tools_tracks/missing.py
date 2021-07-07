from starter2 import *


import three_loopers_1tff as tl

if 'clobber' not in dir():
    clobber=False

frame = 0
if 'hullz' not in dir():
    hullz = {}
    for loop in [tl.looper1, tl.looper2, tl.looper3]:
        hull = CHT.hull_tool(loop)
        hull.make_hulls(frames=[frame])
        hullz[ loop.out_prefix ] = hull

def find_missing_particles(loop,core_id, frame=0):
    nframe=frame
    ms = trackage.mini_scrubber( loop.tr, core_id)
    points = np.column_stack([ ms.this_x[:,nframe], ms.this_y[:,nframe], ms.this_z[:,nframe]])
    final_point = np.column_stack([ ms.mean_xc[-1], ms.mean_yc[-1], ms.mean_zc[-1]])

    #get a rectangel from the data
    left = copy.copy(points.min(axis=0))
    right =copy.copy(points.max(axis=0))
    center = 0.5*(left+right)
    ds = loop.load(frame)
    reg = ds.region(center,left,right)
    other_points = reg['particle_position'].v
    other_ind = reg['particle_index']

    #period shift of particles on disk.  
    #This is not the best method, revise.
    #Shift based on distance of the particle to it's final position:
        #if its more than half the box away from its final place, shift
    for dim in range(3):
        delta = final_point[:,dim]- other_points[:,dim] 
        to_shift = np.abs(delta) > 0.5
        shift = np.zeros_like(delta)
        shift[to_shift] = np.sign( delta[to_shift])

        other_points[to_shift,dim] += shift[to_shift]

    this_hull = hullz[sim_name].hulls[core_id]

    all_points = hullz[sim_name].points_3d[core_id]

    mask_in = CHT.in_hull( other_points, all_points)
    hull_indices = other_ind[mask_in]
    return mask_in, hull_indices


this_simname = 'u202'
loop = tl.looper2   
mask,Core_258_indices = find_missing_particles( loop , 258)
core_list = [2580]
frame_list = loop.tr.frames
derived=[]
color_dict={2580:'r'}
if 'missing_loop' not in dir() :
    #for making
    missing_loop = looper.core_looper(directory= dl.sims[this_simname],
                                     sim_name = this_simname + "_m258",
                                     out_prefix = this_simname,
                                     target_frame = dl.target_frames[this_simname],
                                     frame_list = frame_list,
                                     core_list =  core_list,
                                     fields_from_grid=fields,
                                     derived = derived,
                                     do_shift=False
                                  )
    missing_loop.plot_directory = "/home/dccollins/PigPen"
    missing_loop.target_indices[2580] = Core_258_indices
    #                                 bad_particle_list=dl.bad_particles.get(this_simname,None))
    missing_loop.get_tracks()
    loop_tools.re_shift_snaps( missing_loop)
loop_apps.core_proj_multiple( missing_loop, axis_list=[0], color_dict=color_dict)
