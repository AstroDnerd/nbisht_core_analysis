from starter2 import *
import pyximport; pyximport.install()
import time

def get_current_mask(self):
    """get the particle mask that relates particles for this core_id and this frame
    to the particles in the target_indices from the target_frame"""
    these_pids=self.target_indices.astype('int64')
    if type(these_pids) == yt.units.yt_array:
        these_pids = these_pids.v
    data_region = self.get_region(self.frame)
    mask_to_get=np.zeros(these_pids.shape,dtype='int32')
    print('a1')
    my_indices = sorted(data_region['particle_index'].astype('int64'))
    print('a2')

    found_any, mask = particle_ops.mask_particles_sorted(these_pids,my_indices,mask_to_get)
    print('a3')
    self.mask = mask
    return found_any, mask
def test45():
    this_simname = 'u201'

    output_name = '%s_first_last_t1_nXXX0.h5'%(this_simname)
    new_loop = looper.core_looper(directory= dl.sims[this_simname])
    new_loop.load_loop(output_name)

    ds = new_loop.load(0)

    data_region = ds.all_data() #ds.region([0.5]*3,[0.3]*3,[0.7]*3)
    #data_region = ds.all_data()
    size = 0.25
    left = 0.0
    center = left+size
    right = left+2*size
    center = [center]*3
    left = [left]*3
    right = [right]*3
    #test_region = ds.region(center, left, right)
    test_region = ds.all_data()

    search_for = np.sort(test_region['particle_index'].astype('int64'))[::]
    search_in_unsorted = data_region['particle_index'].astype('int64')[::]
    search_in = np.ascontiguousarray(np.sort(search_in_unsorted))

    mask_to_get = np.zeros_like(search_for, dtype='int32')
    t0 = time.time()
    #found_any, mask1 = particle_ops.mask_particles_sorted_t4(search_for,search_in,mask_to_get)
    t1 = time.time()
    mask_to_get = np.zeros_like(search_for, dtype='int32')
    found_any, mask2 = particle_ops.mask_particles_sorted_t7(search_for,search_in,mask_to_get)
    t2 = time.time()
    s1 = set(search_for.v)
    s2 = set(search_in)
    s3 = s1.intersection(s2)
    t3 = time.time()

    #print("T4 %0.2e %d"%(t1-t0, mask1.sum()))
    print("T5 %0.2e %d"%(t2-t1, mask2.sum()))
    print("Interectin %0.2e"%(t3-t2))


def do_test():
    this_simname = 'u201'

    output_name = '%s_first_last_t1_nXXX0.h5'%(this_simname)
    new_loop = looper.core_looper(directory= dl.sims[this_simname])
    new_loop.load_loop(output_name)

    ds = new_loop.load(0)

    data_region = ds.all_data() #ds.region([0.5]*3,[0.3]*3,[0.7]*3)
    #data_region = ds.all_data()
    size = 0.25
    left = 0.0
    center = left+size
    right = left+2*size
    center = [center]*3
    left = [left]*3
    right = [right]*3
    #test_region = ds.region(center, left, right)
    test_region = ds.all_data()

    search_for = np.sort(test_region['particle_index'].astype('int64'))[::2]
    search_in_unsorted = data_region['particle_index'].astype('int64')[::2]
    search_in = np.ascontiguousarray(np.sort(search_in_unsorted))

    mask_to_get = np.zeros_like(search_for, dtype='int32')
    t0 = time.time()
    found_any, mask1 = particle_ops.mask_particles_sorted_t4(search_for,search_in,mask_to_get)
    t1 = time.time()
    
    print(t0)
    print(t1)
    print("sorted ",t1-t0, " got", mask1.sum())

    mask_to_get = np.zeros_like(search_for, dtype='int32')
    t3 = time.time()
    unsorted = search_in_unsorted
    found_any, mask2 = particle_ops.mask_particles(search_for,unsorted,mask_to_get)
    t4 = time.time()

#   print("unsort ",t4-t3, " got", mask2.sum())
import cProfile
if 1:
    cProfile.run("test45()", "profile_output.txt")
if 1:
    import pstats
    from pstats import SortKey
    p = pstats.Stats('profile_output.txt')
    p.sort_stats(SortKey.TIME)
    p.print_stats(10)

