
from starter2 import *
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
reload(dl)
plt.close('all')
color={'u05':'r','u10':'g','u11':'b'}
color.update({'u201':'r','u202':'g','u203':'b'})

import convex_hull_tools as CHT
reload(CHT)

#import three_loopers_1tff as tl
#core_id = 206
#core_id = 84
#core_id = 112
#core_id = 191
core_id = 84
if 'loop' not in dir():
    file_list=glob.glob("/data/cb1/Projects/P19_CoreSimulations/CoreSets/u203/*c%04d*"%core_id)
    #file_list=glob.glob("/data/cb1/Projects/P19_CoreSimulations/CoreSets/u203/*")
    loop=looper.core_looper(directory=dl.sims['u203'])
    for nfile,fname in enumerate(file_list):
        loop.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    loop.out_prefix='u203_tmp'
    loop.tr.sort_time()


for frame in [0]:
    if 0:
        import three_loopers_1tff as tl
        #
        # Compute overlaps
        #
        ht3 = CHT.hull_tool(tl.looper3)
        for htool in [ht3]: #[ht1, ht2, ht3]:
            htool.overlaps=defaultdict(list)
            htool.make_hulls(frames=[frame])
            for core_1 in htool.cores_used:
                print("overlap li,", core_1)
                for core_2 in htool.cores_used:
                    result = htool.check_hull_overlap(core_1,core_2)
                    if core_1 == core_2:
                        result = -result
                    htool.overlaps[core_1].append(result)
        fractions,cores_84=CHT.get_overlapping_cores(ht3,core_id)
    #cores_84 = np.array([188, 185, 191, 190, 189, 134, 127, 148,  90, 173, 132, 113,  89,
    #                  119, 171, 118, 207,  84, 121, 124, 206, 116, 117, 112, 161, 156,
    #                  178, 166, 110, 163, 182, 157, 105, 164, 109, 106, 162, 158, 165,
    #                  108, 160, 107, 159])



#84, 173, 210
if 1:
    for frame in [0]:# list(range(0,110,10)) + [107]:
        if 0:
            #
            # Compute hulls
            #
            for htool in [ht3]: #[ht1, ht2, ht3]:
                htool.overlaps=defaultdict(list)
                htool.make_hulls(frames=[frame])
        if 1:
            CHT.plot_2d(ht3,core_list = cores_84, accumulate=True,frames=[frame], all_plots=False)
if 0:
    #
    # time series for one core
    #
    for frame in range(0,110,10):
        ht3b = CHT.hull_tool(loop)
        ht3b.make_hulls( core_list = [core_id])
        CHT.plot_2d(ht3b,core_list = [core_id], accumulate=False,frames=[frame], all_plots=True)



if 0:
    reload(trackage)
    #potentially very sweet.
    fig,ax=plt.subplots(1,2)
    thtr = loop.tr
    ms = trackage.mini_scrubber(thtr,core_id)
    if 'bads' not in dir():
        bads = []
    for nnn,p in enumerate(ms.raw_x):
        deelta =  np.abs(p[1:]-p[:-1])
        if np.max(deelta) > 0.5:
            #print("MOOO", np.where( deelta > 0.5))
            ax[0].plot(p,c='k')
            shifted = trackage.shift_6(p)
            #if shifted.min() < -0.8:
            #    bads.append(nnn)
            ax[1].plot(shifted)
            ax[0].plot(shifted)
    fig.savefig('plots_to_sort/shift_debug_c%04d.png'%core_id)

if 0:

    fig,ax=plt.subplots(1,1)
    thtr = loop.tr
    diffz = []
    jumpers = []
    for core_id in np.unique(thtr.core_ids):
        ms = trackage.mini_scrubber(thtr,core_id)
        for arr in [ms.raw_x, ms.raw_y, ms.raw_z]:
            diff = np.abs( arr[:, 1:] - arr[:, :-1])
            diffz.append(diff.max())
            if diff.max() > 0.5:
                jumpers.append(core_id)

    plt.clf()
    plt.hist(diffz)
    plt.savefig('plots_to_sort/max_diff.png')

        #for p in ms.raw_x:
        #    #ax.plot(p)
        #fig.savefig('plots_to_sort/shift_debug_c%04d.png'%core_id)
