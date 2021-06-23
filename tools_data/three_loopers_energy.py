from starter2 import *
reload(looper)
reload(dl)
debug_slice = slice(None)
if 'looper1' not in dir():
    file_list=glob.glob(dl.energy_vorticity_files['u05'])[debug_slice]
    looper1=looper.core_looper(directory=dl.sims['u05'])
    for nfile,fname in enumerate(file_list):
        looper1.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    thtr = looper1.tr
    looper1.core_list = looper1.target_indices.keys()
    looper1.out_prefix='u05'
    thtr.sort_time()

if 'looper2' not in dir():
    looper2=looper.core_looper(directory=dl.sims['u10'])
    file_list=glob.glob(dl.energy_vorticity_files['u10'])[debug_slice]
    for nfile,fname in enumerate(file_list):
        looper2.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    looper2.out_prefix='u10'
    looper2.core_list = looper2.target_indices.keys()

if 'looper3' not in dir():
    looper3=looper.core_looper(directory=dl.sims['u11'])
    file_list=glob.glob(dl.energy_vorticity_files['u11'])[debug_slice]
    for nfile,fname in enumerate(file_list):
        looper3.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    looper3.out_prefix='u11'
    looper3.core_list = looper3.target_indices.keys()
