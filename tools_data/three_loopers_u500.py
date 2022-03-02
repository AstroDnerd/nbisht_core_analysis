
from starter2 import *
reload(looper)
import looper2

if 'loops' not in dir():
    loops={}

for ns,sim_name in enumerate(['u501','u502','u503']):
    print('load',sim_name)
    setname = dl.coresets+"/u500/%s_all_frame_all_prim.h5"%sim_name
    directory = dl.sims[sim_name]
    loops[sim_name] = looper2.load_looper( setname, directory=directory)

