from starter2 import *
from collections import defaultdict
import scipy
import colors

import convex_hull_tools as CHT
import hair_dryer
import otherones
reload(hair_dryer)
import looper2
import three_loopers_six as TL
import otherones
reload(otherones)
import supersets
sim_list=['u602']
if 'big_looper' not in dir():
    big_looper={}
for sim in sim_list:
    if sim in big_looper:
        continue
    which_sim = sim[-1]
    save = {'1': '/data/cb1/Projects/P19_CoreSimulations/CoreSets/a000/a001_all_particles.h5',
                          '2':'/data/cb1/Projects/P19_CoreSimulations/CoreSets/a000/b002_all_particles.h5',
                          '3':'/data/cb1/Projects/P19_CoreSimulations/CoreSets/a000/a003_all_particles.h5'}[which_sim]

    big_looper[sim] = looper2.core_looper2( directory = dl.sims[sim], savefile_only_trackage=save)
    print('make ms, takes about 90 sec')
    ms = trackage.mini_scrubber(big_looper[sim].tr, 0, do_velocity=False)
    big_looper[sim].ms = ms

if 'hull_by_frame' not in dir():
    #
    # Hull volume vs total cell volume
    #
    hull_by_frame = {}
    #looper_list=[MOD.loops['u401']] #,MOD.loops['u402'],MOD.loops['u403']]
    #looper_list=[MOD.loops['u402']] #,MOD.loops['u402'],MOD.loops['u403']]
    #loopers = MOD.loops
#frame_list_1 = [0, 50]
if 1:
    for sim in sim_list:
        if sim in hull_by_frame:
            continue
        loop=TL.loops[sim]
        name = loop.sim_name
        hull_by_frame[name]={}
        if 'frame_list_1' in dir():
            frame_list=frame_list_1
        else:
            frame_list=loop.tr.frames
        #frame_list=[125]
        #frame_list=[0]
        for nframe, frame in enumerate(frame_list):
            hull_by_frame[name][frame]=CHT.hull_tool(loop)
            hull_by_frame[name][frame].make_hulls(frames=[frame])
            hull_by_frame[name][frame].make_overlaps()

if 0:
    import convex_hull_plot2d as CHT2d
    for sim in sim_list:
        for frame in hull_by_frame[sim]:
            fig,ax=plt.subplots(3,1)
            CHT2d.plot_2d(hull_by_frame[sim][frame],core_list=[74],frames=[frame], axis_to_plot=[-1], external_axis=ax)
            fig.savefig('plots_to_sort/others_n%04d.png'%frame)

if 'superset_by_frame' not in dir():
    superset_by_frame={}

if 1:
    for sim in sim_list:
        if sim not in superset_by_frame:
            superset_by_frame[sim] = {}
        for frame in hull_by_frame[sim]:
            if frame in superset_by_frame[sim]:
                continue
            superset_by_frame[sim][frame] = supersets.superset( TL.loops[sim], hull_by_frame[sim][frame])
            superset_by_frame[sim][frame].find()

if 'others' not in dir():
    others={}
reload(otherones)
for sim in sim_list:
    if sim in others:
        continue
    others[sim]={}
    frame_list = sorted( hull_by_frame[sim].keys())
    frame_list = frame_list[-1:]
    core_list=[74]
    for frame in frame_list:
        new_looper = otherones.find_other_ones(sim, hull_by_frame[sim][frame], big_loop=big_looper[sim],frame=frame,
                                               core_list=core_list,superset=superset_by_frame[sim][frame], add_sphere=True)
        others[sim][frame]=new_looper

import otherones_hair
reload(otherones_hair)
if 1:
    for sim in sim_list:
        
        frame_list = sorted( hull_by_frame[sim].keys())
        frame_list=frame_list[-1:]
        fractions=defaultdict(list)
        main_loop =  TL.loops[sim]
        for frame in frame_list:
            other_loop = others[sim][frame]
            if other_loop is None:
                continue
            #if other_loop is not None:
            #    IM = otherones_hair.images(other_loop,main_loop)
            #    IM.run(frames=[0,frame],core_list=core_list,output_prefix='others_n%04d'%frame)
            import xyz_tracks
            reload(xyz_tracks)
            import xyz_twoloop
            reload(xyz_twoloop)
            tracks = xyz_twoloop.tracks(main_loop,other_loop)
            tracks.run(core_list=core_list, output_suffix='others_n%04d'%frame)
if 0:
    for sim in sim_list:
        frame_list = sorted( hull_by_frame[sim].keys())
        fractions=defaultdict(list)
        main_loop =  TL.loops[sim]
        for frame in frame_list:
            other_loop = others[sim][frame]
            if other_loop is None:
                fractions[core_id].append(0)
            else:
                core_list = np.unique(others[sim][frame].tr.core_ids)
                for core_id in core_list:
                    main_pids = main_loop.tr.particle_ids[ main_loop.tr.core_ids == core_id]
                    other_pids = (other_loop.tr.particle_ids[ other_loop.tr.core_ids == core_id])
                    fractions[core_id].append( other_pids.size/main_pids.size)
        fig,ax=plt.subplots(1,1)
        for core_id in fractions:
            plt.plot(fractions[core_id])
        fig.savefig('plots_to_sort/other_fractions_%s.png'%sim)
        print('fart')
                


