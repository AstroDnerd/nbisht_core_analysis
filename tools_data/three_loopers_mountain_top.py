from starter2 import *
reload(looper)
reload(dl)

if 'loops' not in dir():
    loops={}

def load_looper(simname, clobber=False):
    if simname in loops:
        return
    directory = dl.sims[simname]
    savefile =  dl.mountain_top[simname]
    print("Loading loop", simname, savefile)
    loop = looper.core_looper(directory= directory,savefile_only_trackage=savefile)
    loops[simname] = loop

if 'u301' not in loops:
    print('Load u301')
    load_looper('u301')

if 'u302' not in loops:
    print('Load u302')
    load_looper('u302')

if 'u303' not in loops:
    print('Load u303')
    load_looper('u303')

