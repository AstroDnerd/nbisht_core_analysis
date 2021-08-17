"""
This just pulls partilce data from enzo and stores it.
Changing
core_list 
frame_list
fields
changes what gets extracted.
"""
from starter2 import *
import xtra_energy
import data_locations as dl
reload(dl)
#
# set sim
#
if 'this_simname' not in dir():
    this_simname = 'u11'

def three_way_bean():
    from matplotlib.ticker import NullFormatter
    nullfmt = NullFormatter()         # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(1, figsize=(8, 8))

    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    return axScatter,axHistx, axHisty

read_from_disk=True
for this_simname in ['u201','u202','u203']:
    output_name = '%s_first_last_t3_nXXX0.h5'%(this_simname)
    if 'this_looper' not in dir() and not read_from_disk:
        #
        # make first_last object
        #
        directory = dl.sims[this_simname]
        core_list =  [0]
        target = dl.target_frames[this_simname]
        frame_list = [0, target]
        frame_list = list(range(0,target,10)) + [target]

        #frame_list = list(range(0, dl.target_frames[this_simname],10))+[dl.target_frames[this_simname]]
        fields = ['x','y','z','density','velocity_magnitude','magnetic_field_strength', 'velocity_divergence']
        fields += ['velocity_x','velocity_y','velocity_z']
        fields += ['magnetic_field_%s'%s for s in 'xyz']
        #fields += ['PotentialField']
        output_base = "primitive_test"
        derived=[]
        pdb.set_trace()


        this_looper = looper.core_looper(directory= directory,
                                         sim_name = this_simname,
                                         out_prefix = this_simname,
                                         target_frame = dl.target_frames[this_simname],
                                         frame_list = frame_list,
                                         core_list =   core_list,
                                         fields_from_grid=fields,
                                         derived = derived
                                      )
        ds = this_looper.load(frame_list[0])
        ad = ds.all_data() #ds.region([0.5]*3,[0.4]*3,[0.6]*3)
        this_looper.target_indices[0]=ad['particle_index']
        this_looper.get_tracks()
        this_looper.save(output_name)

if 1:

    if 'new_loop_list' not in dir():
        new_loop_list=[]
        for this_simname in ['u201','u202','u203']:
            if this_simname in ['u201']:
                continue
            print("READ FROM DISK")
            output_name = '%s_first_last_t3_nXXX0.h5'%(this_simname)
            new_loop = looper.core_looper(directory= dl.sims[this_simname])
            new_loop.load_loop(output_name)
            new_loop_list.append(new_loop)

def test_rho_c(rho_c, density):


    mask = np.zeros_like(density, dtype='bool')

    for n in range(density.shape[1]):
        mask[:,n] = ( density[:,n] > rho_c ) + mask[:,n-1]

    any_come_back= (density[mask] <= rho_c).any()

    return any_come_back, mask

fig,axlist=plt.subplots(3,1)
ax,ax1,ax2 = axlist
for loop in new_loop_list:
    density = copy.copy(new_loop.tr.track_dict['density'])
    x = copy.copy(new_loop.tr.track_dict['x'])
    y = copy.copy(new_loop.tr.track_dict['y'])
    z = copy.copy(new_loop.tr.track_dict['z'])

    rhos = np.logspace(0,7,30)
    #rhos = [6.5e3]
    nback = np.zeros_like(rhos)
    ncros = np.zeros_like(rhos)
    things=[]
    for nrho,rho_crit in enumerate(rhos):
        I,M = test_rho_c( rho_crit, density)
        things.append([I,M])

        return_particles = np.where( (density < rho_crit) * M )
        niq = np.unique( return_particles[0])
        nback[nrho] = niq.size
        ncros[nrho] = np.unique( np.where( M)[0]).size
        print(rho_crit,niq.size)

    ax.plot( rhos,nback)
    ax1.plot( rhos, ncros)
    ax2.plot( rhos, nback/ncros)
axbonk(ax,xlabel='rho',ylabel='Nback',xscale='log',yscale='log')
axbonk(ax1,xlabel='rho',ylabel='Ncross',xscale='log',yscale='log')
axbonk(ax2,xlabel='rho',ylabel='Ncross',xscale='log',yscale='log')
fig.savefig('plots_to_sort/nback.png')

if 0:
    reload(loop_apps)
    import annotate_particles_3
    reload(annotate_particles_3)
    new_loop.plot_directory = "./plots_to_sort"
    rho_crit=1e5
    particles_to_take = np.max( (density < rho_crit)*M, axis=1)
    print(particles_to_take)
    frames = new_loop.tr.frames
    pos_override = {}
    for nf, frame in enumerate(frames):
        pos_override[frame] = np.column_stack( [ x[particles_to_take,nf], y[particles_to_take,nf], z[particles_to_take,nf]])

    loop_apps.core_proj_multiple(new_loop,axis_list=[2], field='density', frame_list=frames[-1:],
                                                              particles=False, lic=False, fields=False, velocity=False, pos_override=pos_override)
