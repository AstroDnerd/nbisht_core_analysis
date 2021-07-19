from starter2 import *
reload(looper)
reload(dl)
loops={}

def load_looper(simname, clobber=False):
    if simname in loops:
        return
    directory = dl.sims[simname]
    savefile =  dl.mountain_top[simname]
    print("Loading loop", simname, savefile)
    loop = looper.core_looper(directory= directory,savefile=savefile)
    loops[simname] = loop

print('hey')
if 'u301' not in loops:
    load_looper('u301')

