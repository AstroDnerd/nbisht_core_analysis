from starter2 import *
import data_locations as dl
from collections import defaultdict

import davetools
reload(davetools)
import mass_tools
reload(mass_tools)

plt.close('all')

import cProfile
if 0:
    #read test
    import three_loopers_mountain_top as TLM
    outname = "load_profile.txt"
    cProfile.run("TLM.load_looper('u302')", outname)
    #tools/gprof2dot.py -f pstats load_profile.txt | dot -Tpng -o output.png

if 0:
    import tracks_read_write
    reload(tracks_read_write)
    #loop = TLM.loops['u302']
    loop = new_looper
    tracks_read_write.save_loop_only_trackage( loop, 'TRACKAGE_TEST.h5')

if 0:
    import tracks_read_write
    reload(tracks_read_write)

    test_file = 'TRACKAGE_TEST_CORE14.h5'
    savefile = 'TRACKAGE_TEST.h5'
    savefile = 'u302_tracks_only_c0-20.h5'
    savefile = 'u302_zero_test_tracks_only.h5'
    savefile = "u302_zero_test_tracks_only.h5"
    #ktracks_read_write.load_loop_only_trackage( loop, )
    loop3 = looper.core_looper(directory= dl.sims['u302'],savefile_only_trackage=savefile)

if 1:

    mt1 = mass_tools.mass_tool(new_looper)
    mt1.run()

if 1:
    reload(mass_tools)
    fig,ax=plt.subplots(1,1)
    mass_tools.plot_mass_tracks(mt1,ax)
    ax.plot( [0,1],[0,1])
    fig.savefig('plots_to_sort/%s_mass_test_read_test.png'%"U302_c0014")

    plt.close(fig)
    

