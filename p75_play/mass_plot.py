
from starter2 import *
import track_loader as TL
import ucore
#reload(ucore)


import radial_sphere 
reload(radial_sphere)
import radial_plots
reload(radial_plots)


sim_list=['u502']
import track_loader as TL
TL.load_tracks(sim_list)
import monster
monster.load(sim_list)

if 'things' not in dir():
    things={}

for sim in sim_list:
    last_etrack = ucore.etrack_list[-1]
    mon = monster.closet[sim]
    for core_id in mon.this_looper.core_by_mode['A']:
        if core_id in things:
            continue
        prof = radial_sphere.profiler(mon,core_id)
        prof.run(frames='short')
        things[core_id]=prof
        print('did')
    for core_id in things:
        radial_plots.mass_particles(core_id,things[core_id],mon, outname = "%s/bvirial_%s_c%04d"%(plot_dir,sim,core_id))

