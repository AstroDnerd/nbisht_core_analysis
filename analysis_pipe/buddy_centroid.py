

from starter2 import *

import buddy_hair
reload(buddy_hair)
import find_other_cores
reload(find_other_cores)

import track_loader as TL

import json

def buddy_centroid(trackname, thresh=0.05, sink_dict_loc = None):
    track = track_info.tracks[trackname]
    TL.load_tracks([trackname])
    this_looper = TL.tracks[trackname]
    thtr=this_looper.tr
    main_core_list = nar(this_looper.core_list)
    main_core_list.sort()
    mini_scrubbers={}
    for core_id in main_core_list:
        print('buddy_centroid',core_id)
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
        ms.particle_pos(core_id)

        mini_scrubbers[core_id]=ms
    other_cores, shift = find_other_cores.get_other_cores( this_looper, main_core_list, mini_scrubbers, thresh=thresh)
    #get cores with total mass above the mean total mass of all cores
    all_masses = []
    for key_core in main_core_list:
        all_masses.append(mini_scrubbers[key_core].mass_total[-1])
    mean_mass = sum(all_masses)/len(all_masses)
    new_core_list = []
    for key_id in range(len(main_core_list)):
        if all_masses[key_id]>mean_mass:
            new_core_list.append(main_core_list[key_id])
    print("New Core List length: ",len(new_core_list))
    for main_core in new_core_list:
        oc = other_cores[main_core]
        #oc = [375,373,268]
        color_dict = colors.make_core_cmap(oc, cmap = 'winter', seed = -1)
        color_dict[main_core]=[1.0,0.0,0.0]
        core_list = sorted(oc+[main_core])
        suffix = "c%04d"%(main_core)
        what_to_plot = 'centroid'
        what_to_plot = 'xyzt'

        sink_core_dict = None
        if sink_dict_loc !=None:
            with open(sink_dict_loc, 'r') as fp:
                sink_core_dict = json.load(fp)


        buddy_hair.buddy_hair(TL.loops[trackname], core_list=core_list,color_dict=color_dict, what_to_plot=what_to_plot, suffix=suffix, shifter=shift.get(main_core,{}), sink_dict = sink_core_dict)
    
    return new_core_list

