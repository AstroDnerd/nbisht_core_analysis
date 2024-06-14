from starter2 import *
import xtra_energy
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import pcolormesh_helper as pch
import colors
import movie_frames 
G = colors.G
plt.close('all')
import scipy.stats

import radial_sphere 
reload(radial_sphere)
import radial_plots
reload(radial_plots)


sim_list=['u502']
import track_loader as TL
TL.load_tracks(sim_list)
import monster
monster.load(sim_list)

if 'clobber' not in dir():
    clobber=False
if 'thingsp' not in dir() or clobber:
    thingsp={}
    thingss={}
for sim in sim_list:
    core_list = [114]
    mon = monster.closet[sim]
    for core_id in core_list:
        print('ok')
        if core_id in thingsp:
            continue
        prof = radial_sphere.profiler(mon,core_id)
        prof.run()
        thingsp[core_id]=prof
        surf = radial_sphere.surface(mon,core_id)
        surf.run()
        thingss[core_id]=surf

    for core_id in core_list:
        radial_plots.mass_moment(thingsp[core_id], outname = "%s/virial_%s_c%04d"%(plot_dir,sim,core_id))
        #radial_plots.energy(thingsp[core_id], outname = "%s/virial_energy_%s_c%04d"%(plot_dir,sim,core_id))
        #radial_plots.four_way(thingsp[core_id], outname = "%s/virial_quartet_%s_c%04d"%(plot_dir,sim,core_id))
        radial_plots.surface(thingss[core_id], thingsp[core_id],outname = "%s/virial_surface_%s_c%04d"%(plot_dir,sim,core_id))

    




