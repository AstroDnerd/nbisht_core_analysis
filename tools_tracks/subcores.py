from starter2 import *



def thing(this_looper,core_id):
    x = this_looper.tr.c([core_id],'x')
    y = this_looper.tr.c([core_id],'y')
    z = this_looper.tr.c([core_id],'z')
    nx = 128
    ii = (x*nx).astype('int')
    jj = (y*nx).astype('int')
    kk = (z*nx).astype('int')
    i2048= ii +nx*jj+nx*nx*kk
    return i2048



sim_list=['u503']
import three_loopers_u500 as TL

for sim in sim_list:
    i2048=thing(TL.loops[sim],core_id=[9])
    #plt.clf()
    #b = np.sort(i2048[:,0]).flatten()
    #plt.plot( b)
    #plt.savefig('plots_to_sort/dumb1.png')
    counter=5
    for n in range(0,i2048.shape[1],2):
        unique_ids, counts = np.unique( i2048[:,n], return_counts=True)
        n_repeats, n_of_those = np.unique( counts, return_counts=True)
        print( max(n_repeats))

