
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

class box_of_rain():
    def __init__(self):
        self.cores_used=[]
        self.mass_total_hull=[]
        self.mass_other_ones=[]
    def write(self,fname):
        fptr=h5py.File(fname,'w')
        fptr['cores_used']=self.cores_used
        fptr['mass_total_hull']=self.mass_total_hull
        fptr['mass_other_ones']=self.mass_other_ones
        fptr.close()
    def read(self,fname):
        fptr=h5py.File(fname,'r')
        self.cores_used      =fptr['cores_used'][()]     
        self.mass_total_hull =fptr['mass_total_hull'][()]
        self.mass_other_ones =fptr['mass_other_ones'][()]
        fptr.close()

def find_other_ones(new_name, hull_tool,core_list=None,frame=0, superset=None):
    this_loop=hull_tool.this_looper
    big_loop = this_loop.big_loop
    oms = big_loop.ms
    big_density = big_loop.tr.c([0],'density')[:,frame]
    big_cell_volume = big_loop.tr.c([0],'cell_volume')[:,frame]

    if core_list is None:
        core_list = this_loop.core_list

    new_looper = looper2.core_looper2(directory= this_loop.directory,
                             sim_name = new_name,
                             out_prefix = new_name,
                             target_frame = this_loop.target_frame,
                             frame_list = big_loop.frame_list,
                             core_list =  core_list,
                             fields_from_grid=big_loop.fields_from_grid,
                             derived = None, #KLUDGE this is going to be a problem later I fear.
                             do_shift=False
                            )

    box = box_of_rain()
    for core_id in core_list:
        box.cores_used.append(core_id)
        print("Find Otherones c%04d"%core_id)
        ms = trackage.mini_scrubber(this_loop.tr, core_id, do_velocity=False)
        print('   ms done')
        this_hull_points = hull_tool.points_3d[core_id]
        all_particles = big_loop.target_indices
        all_points = np.column_stack( [oms.this_x[:,frame], oms.this_y[:,frame], oms.this_z[:,frame]])

        #cut out the region of interest.  Doesnt help that much.
        this_min = this_hull_points.min(axis=0)
        this_max = this_hull_points.max(axis=0)
        #print('kludge')
        #this_min/=4
        #this_max/=4
        ok = ((all_points >= this_min)*(all_points <= this_max)).all(axis=1)

        #do the work.  Wrap it with a timer
        t0=time.time()

        print('   do hull on %d particles'%ok.sum())
        mask_ok = CHT.in_hull( all_points[ok], this_hull_points)

        #the first mask, ok , cuts out particle that are unrelated.
        #the second mask, mask_ok, is the selection we want out of this subset.
        #mask[mask]*=mask_ok updates the first mask with the second mask.
        mask = copy.copy(ok)
        mask[mask] *= mask_ok

        found_particles = all_particles[mask]

        box.mass_total_hull.append((big_density[mask]*big_cell_volume[mask]).sum())

        t1=time.time()
        print("   time ", t1-t0)

        
        #Remove core particles from found particles.  
        #Core particles defined the hull.  found particles are everything.
        #We want to know about the OtherOnes.
        #PA = all particles
        #P1 = PA[mask] = particles in hull
        #Pc = core particles and neighbors
        #Po = P1 - Pc (Pc contains particles not in Pi, this is fine.)
        #Pcc = P1 intersect Pc (all the core+neighbor particles in the hull)
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

        if 1:
            #Remove Core particles from Otherones
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

        box.mass_other_ones.append((big_density[mask]*big_cell_volume[mask]).sum())
        otherone_particle_ids = all_particles[mask]
        otherone_core_ids = np.ones_like(otherone_particle_ids,dtype='int')*core_id
        otherone_field_values={}
        for field in big_loop.tr.track_dict:
            otherone_field_values[field]=big_loop.tr.track_dict[field][mask,:]
        snap = small_snapshot(particle_ids=otherone_particle_ids,
                              core_ids=otherone_core_ids,
                              frames=big_loop.tr.frames,
                              times = big_loop.tr.times,
                              field_dict=otherone_field_values)
        if new_looper.tr is None:
            new_looper.tr = trackage.track_manager(new_looper)
        new_looper.tr.ingest(snap)
    box.mass_other_ones = nar(box.mass_other_ones)
    box.mass_total_hull = nar(box.mass_total_hull)
    new_looper.box = box
    return new_looper
