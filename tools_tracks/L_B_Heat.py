from starter2 import *

import heat_map
reload(heat_map)


dv = np.zeros([n_cores,n_time])
Beta = np.zeros([n_cores,n_time])
Lmag = np.zeros([n_cores,n_time])
for nc,core_id in enumerate(cores_used):
    Beta[nc,:]=beta_global_array[core_id]
    Lmag[nc,:]=angular_global_array[core_id]
    dv[nc,:]=cell_volume[core_id]
Bbins = np.geomspace( Beta[Beta>0].min(), d.max(),64)
Lbins = np.geomspace( Lmag[Lmag>0].min(), b.max(),64)
hist, xb, yb = np.histogram2d( Beta, Lmag, bins=[Bbins,Lbins],weights=dv)
fig,ax=plt.subplots(1,1)
pch.helper(hist,xb,yb,ax=ax)
fig.savefig('file.png')
