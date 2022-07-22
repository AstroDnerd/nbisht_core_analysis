from starter2 import *



fn = '/data/cb1/Projects/P19_CoreSimulations/new_sims/u26_sphere_amr/ED02_0003/Extra02_0003.cpu0000'
#fn = '/data/cb1/Projects/P19_CoreSimulations/new_sims/u26_sphere_amr/DD0000/data0000.cpu0000'
fptr = h5py.File(fn,'r')
G1 = fptr['Grid00000001']
G2 = fptr['Grid00000002']

fig,ax=plt.subplots(1,2)
ax[0].imshow( G1['PotentialField'][()].sum(axis=1))
ax[1].imshow( G2['PotentialField'][()].sum(axis=1))
plt.savefig('plots_to_sort/grid_direct.png')

