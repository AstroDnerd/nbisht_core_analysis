
if 'set_looper' not in dir():
    savefile='u301_long_pos_only.h5'
    set_looper=looper.core_looper(directory= dl.sims['u301'],savefile_only_trackage=savefile)
    thtr = set_looper.tr
    set_looper.out_prefix='core_13etal'
    thtr.sort_time()

    bad = np.where(thtr.track_dict['density'] <= 0)
    bad_pids = thtr.particle_ids[bad[0]]
    for bad_id in bad_pids:
        print(" strip bad particle ", bad_id)
        n_bad_densities= (thtr.p([bad_id],'density')  <= 0).sum()
        if n_bad_densities == 0:
            print("No bad densities.")
            raise
        particle_index = np.where( thtr.particle_ids == bad_id)
        for field in thtr.track_dict.keys():
            arr = thtr.track_dict[field]
            smaller = np.delete( arr, particle_index,axis=0)
            thtr.track_dict[field]=smaller
        thtr.particle_ids =  np.delete(thtr.particle_ids, particle_index)
        thtr.core_ids =  np.delete(thtr.core_ids, particle_index)




stillbad = np.where(thtr.track_dict['density'] <= 0)
print("STILL BAD", stillbad)

if 0:
    long_tool = lyapunov_tool( set_looper)
    long_tool.run( core_list=[24])#,24,184])


