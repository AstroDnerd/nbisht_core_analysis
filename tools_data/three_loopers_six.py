from starter2 import *
import looper2
reload(looper2)
reload(dl)

if 'loops' not in dir():
    loops={}

for ns,sim_name in enumerate(['u601','u602','u603']):
    print('load',sim_name)
    savefile = dl.coresets['ourset']+'u600/%s_every_10_all_prim.h5'%sim_name
    directory = dl.sims[sim_name]
    loops[sim_name] = looper2.load_looper( savefile, directory=directory)

