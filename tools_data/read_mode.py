
from starter2 import *
from collections import defaultdict

def read(sim_name):

    if 0:
        sim_id = int(sim_name[-1]) #last character is all we need. 501=601
        fname = "browser_data/core_formation_mode_u50%d.h5"%sim_id
        fptr=h5py.File(fname,'r')
        core_ids = fptr['core_ids'][()]
        modes_str = fptr['modes'][()].astype(str)
        modes_tmp = [sss.split(',') for sss in modes_str]
    else:
        sim_id = int(sim_name[-1]) #last character is all we need. 501=601
        fname = "browser_data/Core Browser Plots - u50%s.tsv"%sim_id
        fptr = open(fname)
        lines = fptr.readlines()
        fptr.close()
        core_ids=[]
        modes_str=[]
        for nc,line in enumerate(lines[1:]):
            stuff = line.split('\t')
            core_name = stuff[0]
            core_id = int( core_name.split("c")[1])
            core_ids.append(core_id)
            modes_str.append( stuff[1])
        modes_tmp = [sss.split(',') for sss in modes_str]
        #trim off whitespace. 
        modes_tmp_2 = []
        for core_modes in modes_tmp:
            clean_modes = [mmm.strip() for mmm in core_modes]
            modes_tmp_2.append(clean_modes)
        modes_tmp = modes_tmp_2 

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
            if m.startswith('C'):
                add_cluster=True
            if m.startswith('Alone'):
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


