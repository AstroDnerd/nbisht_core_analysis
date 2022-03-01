from starter2 import *
import looper2
reload(looper)
reload(dl)

if 'loops' not in dir():
    loops={}

if 'u601' not in loops:
    savefile = dl.coresets['ourset']+'u600/u601_every_10_all_prim.h5'
    #savefile = dl.coresets+"/u600/u601_every_10_all_prim.h5"
    print('savefile',savefile)
    loops['u601'] = looper2.load_looper(savefile)
if 'u602' not in loops:
    savefile =dl.coresets['ourset']+ 'u600/u602_every_10_all_prim.h5'
    #savefile =dl.coresets+"/u600/u602_every_10_all_prim.h5"
    loops['u602'] = looper2.load_looper(savefile)
if 'u603' not in loops:
    savefile =dl.coresets['ourset']+ 'u600/u603_every_10_all_prim.h5'
    #savefile =dl.coresets+"/u600/u603_every_10_all_prim.h5"
    loops['u603'] = looper2.load_looper(savefile)
