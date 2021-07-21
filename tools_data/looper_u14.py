
from starter2 import *
reload(looper)
if 'looper14' not in dir():
    this_simname = 'u14'
    file_list=['/data/cb1/Projects/P19_CoreSimulations/CoreSets/MiscTestData/u14_all_primitives_primitives_c0000_nXXXX.h5']
    looper14=looper.core_looper(directory=dl.sims[this_simname])
    for nfile,fname in enumerate(file_list):
        looper14.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
        print( "    ",fname)
    looper14.tr.sort_time()
