from starter2 import *
import xtra_energy
reload(xtra_energy)
import three_loopers_u500 as TL

#ds=TL.loops['u501'].load(111)
ds = yt.load('/data/cb1/Projects/P19_CoreSimulations/new_sims/u26_sphere_amr/DD0001/data0001')
YT_divv = ('gas','divv')

import xtra_operators as xo
def derp(field,data):
    #gi = xtra_energy.grad(data,YT_velocity_x,0)
    gi = xtra_energy.grad(data,YT_velocity_y,1)
    #gi+= xtra_energy.grad(data,YT_velocity_z,2)
    return gi
ds.add_field(YT_divv,derp,validators=[yt.ValidateSpatial(1,[YT_velocity_x, YT_velocity_y, YT_velocity_z])],
             units='1/s',sampling_type='cell')
xtra_energy.add_gravity(ds)
if 0:
    g = ds.index.grids[-1][YT_density_grad_x]
    fig,ax=plt.subplots(1,1)
    ax.imshow(g.sum(axis=1))
    fig.savefig('plots_to_sort/derp.png')

if 1:
    g = ds.index.grids[-1][YT_density_grad_x]
    fig,ax=plt.subplots(1,1)
    ax.imshow(g.sum(axis=1))
    fig.savefig('plots_to_sort/derp1.png')
if 1:
    g = ds.index.grids[0][YT_density_grad_x]
    fig,ax=plt.subplots(1,1)
    ax.imshow(g.sum(axis=1))
    fig.savefig('plots_to_sort/derp0.png')
if 0:
    proj = ds.proj(YT_density_grad_x,1)
    pw=proj.to_pw()
    pw.save('plots_to_sort/herp')
if 0:
    proj = yt.SlicePlot(ds,'x',YT_grav_energy)
    proj.save('plots_to_sort/herp')
    
