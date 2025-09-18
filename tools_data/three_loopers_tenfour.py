from starter2 import *
reload(looper)
reload(dl)
if 'loops' not in dir():
    loops={}

def load_looper(newname, mountain, clobber=False):
    directory = dl.sims[mountain]
    savefile =  dl.mountain_top[mountain]
    print("Loading loop", newname, savefile)
    loop = looper.core_looper(directory= directory,savefile_only_trackage=savefile)
    loops[newname] = loop

if 'u401' not in loops:
    print('Load u401')
    load_looper('u401', 'u301')

if 'u402' not in loops:
    print('Load u402')
    load_looper('u402', 'u302')

if 'u403' not in loops:
    print('Load u403')
    load_looper('u403', 'u303')

for counter in '123':
    mountain_sim = 'u30'+counter
    new_sim = 'u40'+counter
    mountain_top_fname = "datasets_small/%s_mountain_tops_take_9.h5"%mountain_sim
    loops[new_sim].read_targets(mountain_top_fname)
    loops[new_sim].sim_name = new_sim
    loops[new_sim].out_prefix = new_sim

def cut_cores_tenfour(loop):

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



cut={}
for this_simname in ['u401','u402','u403']:
    cut[this_simname]=cut_cores_tenfour( loops[this_simname])





