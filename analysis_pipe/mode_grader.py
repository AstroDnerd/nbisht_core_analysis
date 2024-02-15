from starter2 import *

import find_other_cores
reload(find_other_cores)
import track_loader
import pdb


if 0:
    #
    # Old things that I used to use.
    #
    import anne
    anne.make_hulls()
    if 'mini_scrubbers' not in dir():
        mini_scrubbers={}

    if 1:
        for sim in sim_list:
            if sim in mini_scrubbers:
                continue
            print('mini scrubbers.',sim)
            this_looper=TL.loops[sim]
            mini_scrubbers[sim]={}
            all_cores=np.unique(this_looper.tr.core_ids)
            thtr=this_looper.tr

            for core_id in all_cores:
                ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
                ms.particle_pos(core_id)

                mini_scrubbers[sim][core_id]=ms

    if 'overlaps' not in dir():
        overlaps={}
        ratios={}

    if 1:
        for sim in sim_list:
            ht0=anne.ht[sim]
            loop=TL.loops[sim]
            frame_list=loop.tr.frames
            #frame_list=[0]
            ncores=len(ht0.cores_used)
            nframes=len(frame_list)

            nframes=1
            nframe=0
            if sim not in overlaps:
                overlaps[sim] = np.zeros([ncores,ncores,nframes])
                ratios[sim] = np.zeros([ncores,ncores,nframes])
                htool = anne.ht[sim]
                for nc1, core_id_1 in enumerate(htool.cores_used):
                    for nc2, core_id_2 in enumerate(htool.cores_used):
                        if nc1==nc2:
                            continue
                        a,b=htool.overlaps[core_id_1][nc2], htool.overlaps[core_id_2][nc1] 
                        #please come back to think about this logic.
                        #david, 2022-09-22
                        #2023-01-06 I think it should be zero of either is zero.
                        #The only issue is full containment and zero volume.
                        ratio = max([a,b])
                        rat= sorted( [a,b])
                        if rat[1] == 0:
                            #ratio = max([a,b])
                            ratio = 0
                        else:
                            ratio=rat[0]/rat[1]
                        ratios[sim][nc1,nc2,nframe] = ratio
                        overlaps[sim][nc1,nc2,nframe] = a


def add_cores_to_group(core_id,group=None,buddies=None,all_cores=None):
    group.append(core_id)
    for core_2 in buddies[core_id]:
        if core_2 in all_cores:
            all_cores.remove( core_2)
            add_cores_to_group(core_2, group=group, buddies=buddies,all_cores=all_cores)

import track_loader as TL
def grade_modes(trackname,delta_alone = 0.025):
    track=track_info.tracks[trackname]
    if track.mode_fname == None:
        print("Mode filename not defined.  Please define.  Exiting without doing anything.")
        return 1
    if os.path.exists(track.mode_fname):
        print("Mode file exists, exiting without doing anything.")
        return 0
    TL.load_tracks([trackname])
    this_looper = track_loader.tracks[trackname]
    thtr = this_looper.tr

    temp_mode={}
    buddy_list={}
    mini_scrubbers={}
    for core_id in this_looper.core_list:
        ms =  trackage.mini_scrubber(thtr,core_id, do_velocity=False)
        ms.particle_pos(core_id)
        mini_scrubbers[core_id]  = ms

    temp_mode={'One':[], 'Multi':[]}
    buddy_list={}
    main_core_list=np.unique(this_looper.tr.core_ids)
    min_d,cores_used, distance= find_other_cores.get_all_distances( this_looper, main_core_list, mini_scrubbers)
    distance[ distance==0] = 4000
    dmin = distance.min(axis=2)

    #Temporarily, anything with more than one neighbor is Cluster.
    #After we sort out all the groupins, we'll count.  No sense in repeating machinery.
    for ncore,core_id in enumerate(main_core_list):
        my_args = np.argsort( dmin[ncore])
        my_min = dmin[ncore][my_args]

        near_cores = my_min < delta_alone
        if near_cores.sum() <1:
            my_mode = 'One'
        else:
            my_buddies = cores_used[my_args][near_cores]
            buddy_list[core_id] = my_buddies
            my_mode = 'Multi'
        temp_mode[my_mode].append(core_id)

    A_cores = copy.copy(temp_mode['One'])
    C_cores = set(temp_mode['Multi'])
    groups={}
    group_id=0
    new_modes = []
    while len(C_cores):
        groups[group_id]=[]
        core_id=C_cores.pop()
        add_cores_to_group(core_id,groups[group_id], buddy_list,C_cores)
        group_id+=1

    A_number=0
    B_number=0
    C_number=0
    group_name={}
    for group_id in sorted(groups.keys()):
        if len(groups[group_id])==2:
            group_name[group_id] = "B%d"%B_number
            B_number+=1
        else:
            group_name[group_id] = "C%d"%C_number
            C_number+=1
    for core_id in main_core_list:
        name=None
        for group_id in sorted(groups.keys()):
            if core_id in groups[group_id]:
                name = group_name[group_id]
        if name is None:
            if core_id not in A_cores:
                print("serious error.")
                raise
            name = 'A%d'%A_number
            A_number+=1
        new_modes.append(name)

    print(new_modes)
    fptr=h5py.File( track.mode_fname, 'w')
    try:
        #pdb.set_trace()    
        asciiList = [n.encode("ascii","ignore") for n in new_modes]
        fptr.create_dataset('modes', (len(asciiList),1),'S10',asciiList)
        fptr['core_ids']=main_core_list
    except:
        raise
    finally:
        fptr.close()

