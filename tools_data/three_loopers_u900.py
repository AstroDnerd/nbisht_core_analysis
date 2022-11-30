
from starter2 import *
reload(looper)
import looper2

if 'loops' not in dir():
    loops={}

for ns,sim_name in enumerate(['u902']):
    print('load',sim_name)
    setname = dl.coresets['ourset']+"/u900/%s_vel_grad.h5"%sim_name
    directory = dl.sims[sim_name]
    loops[sim_name] = looper2.load_looper( setname, directory=directory)

