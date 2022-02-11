
from starter2 import *
reload(looper)
import looper2

if 'loops' not in dir():
    loops={}

for sim_name in ['u501','u502','u503']:
    print('load',sim_name)
    loops[sim_name] = looper2.load_looper( dl.u500[sim_name])
