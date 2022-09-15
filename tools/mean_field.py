from starter2 import *

import three_loopers_u500 as TL
sim_list=['u501','u502','u503']
if 'B' not in dir():
    print("Compute Mean Field.")
    B = {}
    for ns,sim in enumerate(sim_list):
        loop = TL.loops[sim]
        ds = loop.load(0)
        ad = ds.all_data()
        Bx = ad['magnetic_field_x'].mean()
        By = ad['magnetic_field_y'].mean()
        Bz = ad['magnetic_field_z'].mean()
        B['u50%d'%(ns+1)] = np.sqrt(Bx*Bx+By*By+Bz*Bz)
        B['u60%d'%(ns+1)] = np.sqrt(Bx*Bx+By*By+Bz*Bz)


    print(B)
