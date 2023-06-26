

from starter2 import *

import buddy_hair
reload(buddy_hair)
import find_other_cores
reload(find_other_cores)

import track_loader as TL


def buddy_centroid(trackname, thresh=0.05):
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
    for main_core in main_core_list:
        oc = other_cores[main_core]
        #oc = [375,373,268]
        color_dict = colors.make_core_cmap(oc, cmap = 'winter', seed = -1)
        color_dict[main_core]=[1.0,0.0,0.0]
        core_list = sorted(oc+[main_core])
        suffix = "c%04d"%(main_core)
        what_to_plot = 'centroid'
        what_to_plot = 'xyzt'

        buddy_hair.buddy_hair(TL.loops[trackname], core_list=core_list,color_dict=color_dict, what_to_plot=what_to_plot, suffix=suffix, shifter=shift.get(main_core,{}))

