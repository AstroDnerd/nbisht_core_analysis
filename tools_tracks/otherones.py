
from starter2 import *
import looper2
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
import convex_hull_tools as CHT
reload(CHT)
import time
reload(trackage)

bucket={}
class small_snapshot():
    def __init__(self, particle_ids=None, core_ids=None, frames=None, field_dict=None,times=None):
        self.ind=particle_ids
        self.core_ids=core_ids
        self.frame=frames
        self.field_values=field_dict
        self.time = times
def find_other_ones(new_name, hull_tool,core_list=None,frame=0, superset=None):
    this_loop=hull_tool.this_looper
    other_loop = this_loop.big_loop
    oms = other_loop.ms

    if core_list is None:
        core_list = this_loop.core_list

    new_looper = looper2.core_looper2(directory= this_loop.directory,
                             sim_name = new_name,
                             out_prefix = new_name,
                             target_frame = this_loop.target_frame,
                             frame_list = other_loop.frame_list,
                             core_list =  core_list,
                             fields_from_grid=other_loop.fields_from_grid,
                             derived = None, #KLUDGE this is going to be a problem later I fear.
                             do_shift=False
                            )

    for core_id in core_list:
        ms = trackage.mini_scrubber(this_loop.tr, core_id, do_velocity=False)
        this_hull_points = hull_tool.points_3d[core_id]
        all_particles = other_loop.target_indices
        all_points = np.column_stack( [oms.this_x[:,frame], oms.this_y[:,frame], oms.this_z[:,frame]])

        #cut out the region of interest.  Doesnt help that much.
        this_min = this_hull_points.min(axis=0)
        this_max = this_hull_points.max(axis=0)
        #print('kludge')
        #this_min/=4
        #this_max/=4
        print(this_min,this_max)
        #ok = ((all_points >= this_min)*(all_points <= this_max)).all(axis=1)


        #do the work.  Wrap it with a timer
        t0=time.time()

        print('do the hull')
        mask = CHT.in_hull( all_points, this_hull_points)
        found_particles = all_particles[mask]

        t1=time.time()
        print("time ", t1-t0)

        
        #Remove core particles from found particles.  
        #Core particles defined the hull.  found particles are everything.
        #We want to know about the OtherOnes.
        #PA = all particles
        #P1 = PA[mask] = particles in hull
        #Pc = core particles and neighbors
        #Po = Pi - Pc (Pc contains particles not in Pi, this is fine.)
        #Pcc = Pi intersect Pc (all the core+neighbor particles in the hull)
        #PB = [Po] + [Pcc] Make one long array
        #BB = [T]  + [F]   A boolean array along side Pb
        #S1 = Sort(P1)
        #R1 = Sort(S1) #return sort for P1
        #mask contains all particles in the hull.  Turn off the core particles
        # Sort(P1) parallels Sort(BB).  
        # Since P1 may not be sorted, sort BB like P1 with R1
        # mask[mask] *= BB[S1][R1] first sort BB, then sort like P1.

        this_particles = this_loop.target_indices[ core_id]
        if superset is not None:
            my_superset=superset.set_by_core[core_id]
            for other_core in superset.supersets[my_superset]:
                this_particles = np.concatenate([this_particles,
                                                 this_loop.target_indices[other_core]])
        #neighborhoods go here
        this_set = set(this_particles)
        found_set = set(found_particles)
        this_only = found_set.intersection(this_set)
        otherones = found_set - this_only
        otherones_mask = np.concatenate([np.ones(len(otherones),dtype='bool'), 
                                         np.zeros(len(this_only),dtype='bool')])
        particles_both = np.concatenate([nar(list(otherones)), nar(list(this_only))])
        sort_both = np.argsort(particles_both)
        sort1 = np.argsort(found_particles)
        sort_back = np.argsort(sort1) 
        mask[ mask] *= otherones_mask[sort_both][sort_back]

        otherone_particle_ids = all_particles[mask]
        otherone_core_ids = np.ones_like(otherone_particle_ids,dtype='int')*core_id
        otherone_field_values={}
        for field in other_loop.tr.track_dict:
            otherone_field_values[field]=other_loop.tr.track_dict[field][mask,:]
        snap = small_snapshot(particle_ids=otherone_particle_ids,
                              core_ids=otherone_core_ids,
                              frames=other_loop.tr.frames,
                              times = other_loop.tr.times,
                              field_dict=otherone_field_values)
        if new_looper.tr is None:
            new_looper.tr = trackage.track_manager(new_looper)
        new_looper.tr.ingest(snap)
    return new_looper
