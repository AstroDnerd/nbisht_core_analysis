from starter2 import *
import data_locations as dl
from collections import defaultdict

import davetools
reload(davetools)
import mass_tools
reload(mass_tools)

plt.close('all')


if 0:
    import tracks_read_write
    reload(tracks_read_write)

    savefile = 'TRACKAGE_TEST.h5'
    loop = looper.core_looper(directory= dl.sims['u302'],savefile_only_trackage=savefile)

    mt1 = mass_tools.mass_tool(loop)
    mt1.run(core_list = [14])

if 0:
    reload(mass_tools)
    fig,ax=plt.subplots(1,1)
    mass_tools.plot_mass_tracks(mt1,ax)
    fig.savefig('plots_to_sort/%s_mass_test_read_test.png'%mt1.this_looper.out_prefix)
    plt.close(fig)

if 1:
    import three_loopers_mountain_top as TLM
    reload(TLM)
if 1:
    mt_dict={}
    for this_simname in ['u301','u302','u303']:

        mt_dict[this_simname] = mass_tools.mass_tool(TLM.loops[this_simname])
        mt_dict[this_simname].run()
if 1:
    for this_simname in ['u301','u302','u303']:
        fig,ax = plt.subplots(1,1)
        mass_tools.plot_mass_tracks(mt_dict[this_simname], ax)
        fig.savefig('plots_to_sort/%s_mass_time.png'%this_simname)
        plt.close(fig)
if 0:
    import three_loopers_1tff as tl
    if 'clobber' not in dir():
        clobber=True
    if 'mass_tool1' not in dir() or clobber:
        mass_tool1=mass_tools.mass_tool(tl.looper1)
        mass_tool1.run()


    if 'mass_tool2' not in dir() or clobber:
        mass_tool2=mass_tools.mass_tool(tl.looper2)
        mass_tool2.run()
    if 'mass_tool3' not in dir() or clobber:
        mass_tool3=mass_tools.mass_tool(tl.looper3)
        mass_tool3.run()

    fig,ax=plt.subplots(1,3, figsize=(12,4))
    axes=ax.flatten()


    for nt,tool in enumerate([mass_tool1,mass_tool2,mass_tool3]):
        mass_tools.plot_mass_tracks(tool, axes[nt])

    axes[0].set_ylabel(r'$M/\bar{M_{\rm{early}}}$')
    fig.colorbar(ploot)
    fig.savefig(outname)
    print(outname)

