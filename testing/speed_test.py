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
import cProfile
reload(looper)
reload(dl)
#
# set sim
#
this_simname = 'u05'
if 0:
    """Large test.  Works pretty nicely now."""
    core_list =  [0] 
    frame_list = [1]#[0]#+list(range(10,130,10))+[125]
    fields = ['x','y','z','density']
    output_base = "speed_test"
    derived=[]


    this_looper = looper.core_looper(directory= dl.sims[this_simname],
                                     sim_name = this_simname,
                                     out_prefix = this_simname,
                                     target_frame = dl.target_frames[this_simname],
                                     frame_list = frame_list,
                                     core_list =  core_list,
                                     fields_from_grid=fields,
                                     derived = derived
                                  )

    import cProfile
    ds = this_looper.load(frame_list[0])
    #ad = ds.all_data() 
    ad = ds.region([0.5]*3,[0.4]*3,[0.6]*3)
    ad = ds.region([0.5]*3,[0.25]*3,[0.75]*3)
    this_looper.target_indices[0]=ad['particle_index']
    fname = "current_mask_0.txt"
    #cProfile.run("test_snapshot.get_current_mask()",fname)
    cProfile.run("this_looper.get_tracks()", fname)
    #visualize this with:
    #gprof2dot.py -f pstats output.pstats | dot -Tpng -o output.png
    # https://stackoverflow.com/questions/843671/profiling-in-python-who-called-the-function
    #this_looper.get_tracks()
    #this_looper.save(output_name)

if 1:
    """Check that we get the same answer again."""
    this_simname = 'u05'
    core_list =  [10] 
    frame_list = [125]#[0]#+list(range(10,130,10))+[125]
    fields = ['x','y','z','density']
    output_base = "core10b"
    derived=[]
    new_c10 = looper.core_looper(directory= dl.sims[this_simname],
                                     sim_name = this_simname,
                                     out_prefix = this_simname,
                                     target_frame = dl.target_frames[this_simname],
                                     frame_list = frame_list,
                                     core_list =  core_list,
                                     fields_from_grid=fields,
                                     derived = derived
                                  )

    new_c10.get_target_indices(h5_name=dl.peak_list[this_simname],
                                     bad_particle_list=dl.bad_particles[this_simname])
    cProfile.run("new_c10.get_tracks()", "profile_core10.profile")

    if 'actual_core_10' not in dir():
        file_list=['/data/cb1/Projects/P19_CoreSimulations/CoreSets/u05_every_ten/all_primitives_c0010.h5']
        looper_c10=looper.core_looper(directory=dl.sims['u05'])
        for nfile,fname in enumerate(file_list):
            looper_c10.load_loop(fname)
        looper_c10.out_prefix='fiducial_c10'

p1=looper_c10.snaps[125][10].pos
p2=new_c10.snaps[125][10].pos
print("Should be zero: particle positions", np.abs(p1-p2).sum())
d1a =looper_c10.tr.c(10,'density')[:,0]
d2a =new_c10.tr.c(10,'density')[:,0]
def drf(field):
    d1 =looper_c10.tr.c(10,field)[:,0]
    d2 =new_c10.tr.c(10,field)[:,0]
    print("Density should be zero",np.abs(d1-d2).sum())
    print("                out of",np.abs(d1+d2).sum())

if 0:
    import pstats
    from pstats import SortKey
    p = pstats.Stats(fname)
    p.sort_stats(SortKey.TIME)
    p.print_stats(10)
