from starter2 import *

import three_loopers_six as TL
sim_list=['u601','u602','u603']

if 'ct' not in dir():
    ct = {}
    for this_simname in TL.loops:

        ct[this_simname] = close_tool.close_tool( TL.loops[this_simname])
        ct[this_simname].make_distance()

#from analysis_block import *
if 1:
    fig,ax=plt.subplots(1,3)
    for ns,sim in enumerate(TL.loops):
        b = ct[sim].distance_matrix + 0
        #b.sort(axis=0)
        args = np.argsort( b[:,1])

        for nr,row in enumerate(b[args]):
            the_row = row[1:4]
            ax[ns].scatter( [nr]*the_row.size, the_row, c=['r','g','b'])
        axbonk(ax[ns],yscale='log')
    fig.savefig('plots_to_sort/mut.png')
if 0:
    fig,ax=plt.subplots(1,3)
    for ns,sim in enumerate(TL.loops):
        b = ct[sim].distance_matrix + 0
        #b.sort(axis=0)
        args = np.argsort( b[:,1])

        for nr,row in enumerate(b[args]):
            ax[ns].scatter( [nr]*row.size, row)
        axbonk(ax[ns],yscale='log')
    fig.savefig('plots_to_sort/wut.png')

if 0:
    fig,ax=plt.subplots(1,1)
    for row in ct['u601'].distance_matrix:
        ax.plot(row)
    fig.savefig('plots_to_sort/wut.png')

if 0:
    fig,ax=plt.subplots(1,3)
    for ns,sim in enumerate(TL.loops):
        b = ct[sim].distance_matrix + 0
        b.sort(axis=0)

        bins = np.geomspace( b[b>0].min(), b.max(),16)
        #g=ax[ns].hist( b[:,0], histtype='step',color='r',bins=bins)
        ax[ns].hist( b[:,1], histtype='step',color='g',bins=bins)
        ax[ns].hist( b[:,2], histtype='step',color='b',bins=bins)
        ax[ns].hist( b[:,3], histtype='step',color='r',bins=bins)
        axbonk(ax[ns],xscale='log')
    fig.savefig('plots_to_sort/friends.png')


