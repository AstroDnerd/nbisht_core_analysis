from starter2 import *


import three_loopers_1tff as tl
reload(loop_apps)

if 0:
    #does not work on the three_looper_1tff loopers
    loop_tools.re_shift_snaps(tl.looper2)

if 0:
    tl.looper2.plot_directory = './plots_to_sort'
    loop_apps.core_proj_multiple( tl.looper2, core_list=[258], axis_list=[0], color_dict={258:'r'}, 
                                 frame_list = tl.looper2.tr.frames,
                                 #frame_list = [118],
                                zoom=True, only_sphere=True)# tl.looper2.tr.frames)


if 'clobber' not in dir():
    clobber=False

frame = 0
import convex_hull_tools as CHT
if 'hullz' not in dir():
    hullz = {}
    for loop in [ tl.looper2]:
        hull = CHT.hull_tool(loop)
        hull.make_hulls(frames=[frame])
        hullz[ loop.out_prefix ] = hull

def find_missing_particles(loop,hull_tool,core_id, frame=0):
    nframe=frame
    ms = trackage.mini_scrubber( loop.tr, core_id, do_velocity=False)
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

    this_hull = hull_tool.hulls[core_id]

    all_points = hull_tool.points_3d[core_id]

    mask_in = CHT.in_hull( other_points, all_points)
    hull_indices = other_ind[mask_in]
    return mask_in, hull_indices


this_simname = 'u202'
loop = tl.looper2   
mask,Core_258_indices = find_missing_particles( loop , hullz[this_simname], 258)
core_list = [2580]
frame_list = loop.tr.frames
derived=[]
fields = ['x','y','z','density']
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

if 1:
    loop_apps.core_proj_multiple( missing_loop, axis_list=[0], color_dict=color_dict)

if 0:
    ds = missing_loop.load( dl.target_frames[this_simname])
    #ms = trackage.mini_scrubber( missing_loop.tr, 2580, do_velocity=False)
    argmax = np.argmax(missing_loop.tr.track_dict['density'][:,-1])
    xmax = missing_loop.tr.track_dict['x'][argmax,-1]
    ymax = missing_loop.tr.track_dict['y'][argmax,-1]
    zmax = missing_loop.tr.track_dict['z'][argmax,-1]
    x    = missing_loop.tr.track_dict['x'][:,-1]
    y    = missing_loop.tr.track_dict['y'][:,-1]
    z    = missing_loop.tr.track_dict['z'][:,-1]

    c = nar([xmax,ymax,zmax])


    pos = np.column_stack( [x,y,z])
    ok =  x > c[0]-1./16
    ok *= x < c[0]+1./16
    ok =  y > c[1]-1./16
    ok *= y < c[1]+1./16
    ok =  z > c[2]-1./16
    ok *= z < c[2]+1./16

    pos = pos[ok]

    SL = yt.SlicePlot( ds,'x','density',center=c)
    SL.annotate_these_particles2(1.0, col='r', positions=pos)
    SL.zoom(16)
    SL.set_cmap('density','Greys')

    loop=tl.looper2
    x2    = loop.tr.c([258],'x')[:,-1]
    y2    = loop.tr.c([258],'y')[:,-1]
    z2    = loop.tr.c([258],'z')[:,-1]

    SL.annotate_these_particles2(1.0, col='g', positions=np.column_stack([x2,y2,z2]), p_size=10)
    SL.set_axes_unit('code_length')


    SL.save('plots_to_sort/tmp')



