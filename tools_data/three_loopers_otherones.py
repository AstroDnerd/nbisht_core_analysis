from starter2 import *
import looper2
reload(looper)
reload(dl)

if 'loops' not in dir():
    loops={}

if 'x001' not in loops:
    savefile = dl.coresets+"x600/otherones_u401_no_neighbor.h5"
    loops['x001'] = looper2.load_looper(savefile)
if 'x002' not in loops:
    savefile = dl.coresets+"x600/otherones_u402_no_neighbor.h5"
    loops['x002'] = looper2.load_looper(savefile)
if 'x003' not in loops:
    savefile =dl.coresets+ "x600/otherones_u403_no_neighbor.h5"
    loops['x003'] = looper2.load_looper(savefile)
