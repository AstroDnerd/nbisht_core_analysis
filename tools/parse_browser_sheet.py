from starter2 import *
import pandas as pd

sim_list=['u501','u502','u503']
for sim in sim_list:
    core_ids = []
    modes = []
    if 0:
        D1 = pd.read_excel("browser_data/core_browser_plots.xlsx", sheet_name=sim)
        for nc,core_name in enumerate(D1['Core']):
            core_id = int( core_name.split("c")[1])
            core_ids.append(core_id)
            modes.append( D1['Mode'][nc])
        print(modes)
    else:
        fname = "browser_data/Core Browser Plots - %s.tsv"%sim
        fptr = open(fname)
        lines = fptr.readlines()
        fptr.close()
        for nc,line in enumerate(lines[1:]):
            stuff = line.split('\t')
            core_name = stuff[0]
            core_id = int( core_name.split("c")[1])
            core_ids.append(core_id)
            modes.append( stuff[1])

    fptr = h5py.File('browser_data/core_formation_mode_%s.h5'%sim,'w')
    try:
        argsort = np.argsort(core_ids)
        fptr['core_ids'] = core_ids
        fptr['modes'] = modes
    except:
        raise
    finally:
        fptr.close()





