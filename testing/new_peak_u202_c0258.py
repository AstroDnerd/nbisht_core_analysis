from starter2 import *

if 'this_simname' not in dir():
    this_simname = 'u302'

def get_density(peak, test_main):
    X = test_main['x'].to('code_length').v
    Y = test_main['y'].to('code_length').v
    Z = test_main['z'].to('code_length').v
    tolerance = 0.5/2048
    peak = ( np.abs(X - peak[0])<tolerance)*( np.abs(Y-peak[1])<tolerance)*(np.abs(Z-peak[2])<tolerance)
    return test_main['density'][peak]

def get_mountain_top(ds,peak, radius=1e-2):
    test_sphere = ds.sphere(peak,radius)
    peak_density = get_density(peak,test_sphere).v
    min_density = peak_density**(3./4)
    mountain_top = Clump(test_sphere,('gas','density'))
    nsteps=10
    step=(peak_density/min_density)**(1./nsteps)
    find_clumps(mountain_top, min_density,peak_density,step)
    #mountain_top.find_children( min_density, peak_density)
    for leaf in mountain_top.leaves:
        min_radius=leaf['radius'].min()
        if min_radius < 1/2048:
            break
    return leaf
def proj_mountain_top(ds,peak, radius=1e-2):
    test_sphere = ds.sphere(peak,radius)
    leaf = get_mountain_top(ds,peak,radius)

    proj = ds.proj('density',0,center=peak,data_source=test_sphere)
    pw = proj.to_pw()
    pw.set_cmap('density','Greys')

    #pw.annotate_clumps([master_clump]+master_clump.leaves)
    #pw.annotate_these_particles2(1.0,col='r',positions= test_sphere['particle_position'])
    pw.annotate_clumps([leaf])
    pw.zoom(0.5/radius)
    pw.set_axes_unit('code_length')
    print(pw.save('plots_to_sort/test.png'))




if 1:
    frame = 118
    ds = yt.load("%s/DD%04d/data%04d"%(dl.sims[this_simname],frame,frame))
    peak_u202_c0258 = np.array([0.89331055, 0.1159668 , 0.4440918 ])
    peak_u202_c0089 = np.array([0.08569336, 0.12866211, 0.20825195])
    this_peak = peak_u202_c0258
    this_peak = peak_u202_c0089
    this_core = 8901
    mountain_top = get_mountain_top(ds,this_peak)
    proj_mountain_top(ds,this_peak)
    mountain_top_indices = mountain_top['particle_index']
    core_list = [this_core]
    frame_list = loop.tr.frames
    derived=[]
    fields = ['x','y','z','density']
    color_dict={this_core:'r'}
if 0:
    if 'top_loop' not in dir() :
        #for making
        pdb.set_trace()
        top_loop = looper.core_looper(directory= dl.sims[this_simname],
                                         sim_name = this_simname + "_top_89",
                                         out_prefix = this_simname,
                                         target_frame = dl.target_frames[this_simname],
                                         frame_list = frame_list,
                                         core_list =  core_list,
                                         fields_from_grid=fields,
                                         derived = derived,
                                         do_shift=False
                                      )
        top_loop.plot_directory = "/home/dccollins/PigPen"
        top_loop.target_indices[this_core] = mountain_top_indices
        #                                 bad_particle_list=dl.bad_particles.get(this_simname,None))
        top_loop.get_tracks()
        loop_tools.re_shift_snaps( top_loop)
    
if 0:
    loop_apps.core_proj_multiple( top_loop, core_list=[this_core], axis_list=[0], color_dict={2581:'r'}, 
                                 #frame_list = tl.looper2.tr.frames,
                                 #frame_list = [118],
                                zoom=True, only_sphere=True)# tl.looper2.tr.frames)
