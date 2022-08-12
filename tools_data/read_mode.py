
from starter2 import *
from collections import defaultdict

def read(sim_name):

    sim_id = int(sim_name[-1]) #last character is all we need. 501=601
    fname = "browser_data/core_formation_mode_u50%d.h5"%sim_id
    fptr=h5py.File(fname,'r')
    core_ids = fptr['core_ids'][()]
    modes_str = fptr['modes'][()].astype(str)
    modes_tmp = [sss.split(',') for sss in modes_str]
    modes=[]
    unique_modes=[]

    core_by_mode=defaultdict(list)

    for nm, mode in enumerate(modes_tmp):
        add_binary=False
        add_cluster=False
        one=False
        shard=False
        merge=False
        for mmm in mode:
            m = mmm.strip()
            if m not in unique_modes:
                unique_modes.append(m)
            if m.startswith('B'):
                add_binary=True
            if m.startswith('S'):
                add_cluster=True
            if m.startswith('One'):
                one=True
            if m.startswith('Shard'):
                shard=True
            if m.startswith('Merger'):
                merge=True
            core_by_mode[m].append( core_ids[nm])
        if add_binary:
            mode.append('Binary')
            core_by_mode['Binary'].append( core_ids[nm])
        if add_cluster:
            mode.append('Cluster')
            core_by_mode['Cluster'].append( core_ids[nm])
        #this is already done.
        #if one:
        #    core_by_mode['One'].append( core_ids[nm])
        #if add_binary + add_cluster + one == False:
        #    pdb.set_trace()
            
        modes.append(mode)
    for F in core_by_mode:
        core_by_mode[F] = nar( core_by_mode[F])

    mode_dict = dict(zip(core_ids,modes))
    return {'core_ids':core_ids, 'modes':modes, 'mode_dict':mode_dict, 'core_by_mode':core_by_mode, 'unique_modes':unique_modes}


