


from starter2 import *
import data_locations as dl
from collections import defaultdict

import davetools
reload(davetools)
import mass_tools
reload(mass_tools)
import looper2


plt.close('all')

sim_list=['u501']
if 'u501' not in these_loops:
    these_loops={}
    print('read u501')
    these_loops['u501']=looper2.load_looper('u501_all_frame_all_prim.h5')

if 'mtd' not in dir():
    mtd={}
for this_simname in sim_list:
    if this_simname not in mtd:
        mtd[this_simname]=mass_tools.mass_tool(these_loops[this_simname])
        mtd[this_simname].run()


for this_simname in sim_list:
    fig,ax=plt.subplots(1,1)
    mass_tools.plot_mass_tracks(mtd[this_simname],ax)
    fig.savefig('plots_to_sort/thing.pdf')



