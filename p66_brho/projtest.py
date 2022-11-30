
from starter2 import *


dirname = '/data/cb1/Projects/P19_CoreSimulations/new_sims/u28_blank/'
frame = 0
ds = yt.load('%s/DD%04d/data%04d'%(dirname,frame,frame))

L = nar([0.25,0.25,0.25])
R = nar([0.75,0.75,0.75])
c=0.5*(L+R)
reg = ds.region(c,L,R)
pproj = yt.ProjectionPlot( ds, 0, ('gas','density'), weight_field=('gas','cell_volume'), data_source=reg)
dproj = ds.proj( ('gas','density'),0, data_source=reg)#, weight_field=('gas','cell_volume'))
frb=proj.frb
print('word',(dproj['density']*dproj['dx']*dproj['dy']).sum())
print('word',(reg['cell_mass'].sum()))


