
reload(CHT)
if 1:
    for this_simname in sim_list:
        #new_super = supersets.superset(  MOD.loops[this_simname], ht_set[this_simname])
        new_super.find()
        htool=ht_set[this_simname]
        htool.this_looper.out_prefix='u401_c0124etal'
        cores_used=nar(htool.cores_used)
        ooo = nar(htool.overlaps[124])
        overlaps = cores_used[ooo > 0]
        overlaps = np.append(overlaps,124)
        color_dict = colors.make_core_cmap(overlaps, cmap = 'tab20c', seed = -1)
        frame_list = MOD.loops[this_simname].frame_list[0:2]
        CHT.plot_2d(htool,core_list=overlaps,
                    frames=frame_list,accumulate=True,label_cores=[0], prefix = "S%02d"%nset,
                   color_dict=color_dict)

import supersets
reload(supersets)
if 0:
    for this_simname in sim_list:
        new_super = supersets.superset(  MOD.loops[this_simname], ht_set[this_simname])
        new_super.find()
        #supersets = st_set[this_simname]
        #supersets.supersets
        for nset,this_superset in enumerate(new_super.supersets):
            print(this_superset)
            #if 124 not in this_superset:
            #    continue
            color_dict = colors.make_core_cmap(this_superset, cmap = 'tab20c', seed = -1)
            frame_list = MOD.loops[this_simname].frame_list[0:1]
            CHT.plot_2d(ht_set[this_simname],core_list=this_superset,
                        frames=frame_list,accumulate=True,label_cores=[0], prefix = "S%02d"%nset,
                       color_dict=color_dict)
