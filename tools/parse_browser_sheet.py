from starter2 import *
import pandas as pd

sim_list=['u501','u502','u503']
for sim in sim_list:
    D1 = pd.read_excel("browser_data/core_browser_plots.xlsx", sheet_name=sim)
    core_ids = []
    modes = []
    for nc,core_name in enumerate(D1['Core']):
        core_id = int( core_name.split("c")[1])
        core_ids.append(core_id)
        modes.append( D1['Mode'][nc])
    print(modes)
    fptr = h5py.File('browser_data/core_formation_mode_%s.h5'%sim,'w')
    try:
        argsort = np.argsort(core_ids)
        fptr['core_ids'] = core_ids
        fptr['modes'] = modes
    except:
        raise
    finally:
        fptr.close()





