

from starter2 import *

import buddy_hair
reload(buddy_hair)
import find_other_cores
reload(find_other_cores)
import close_tool
reload(close_tool)
import three_loopers_u500 as TL
sim_list = ['u501','u502','u503']
#sim_list = ['u502']

if 'mini_scrubbers' not in dir():
    print('mini scrubbers.')
    mini_scrubbers={}
    for sim in ['u501','u502','u503']:
        this_looper=TL.loops[sim]
        mini_scrubbers[sim]={}
        all_cores=np.unique(this_looper.tr.core_ids)
        thtr=this_looper.tr

        for core_id in all_cores:
            print('ms',core_id)
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            ms.particle_pos(core_id)

            mini_scrubbers[sim][core_id]=ms


for sim in sim_list:
    this_looper = TL.loops[sim]
    main_core_list = np.unique(TL.loops[sim].tr.core_ids)
    main_core_list.sort()
    #main_core_list = [191, 74, 112, 368]
    #main_core_list=[368]
    #main_core_list=[74]
    other_cores, shift = find_other_cores.get_other_cores( this_looper, main_core_list, mini_scrubbers[sim])
    for main_core in main_core_list:
        oc = other_cores[main_core]
        #oc = [375,373,268]
        color_dict = colors.make_core_cmap(oc, cmap = 'winter', seed = -1)
        color_dict[main_core]=[1.0,0.0,0.0]
        core_list = sorted(oc+[main_core])
        suffix = "c%04d"%(main_core)
        what_to_plot = 'centroid'
        what_to_plot = 'xyzt'

        buddy_hair.buddy_hair(TL.loops[sim], core_list=core_list,color_dict=color_dict, what_to_plot=what_to_plot, suffix=suffix, shifter=shift.get(main_core,{}))

