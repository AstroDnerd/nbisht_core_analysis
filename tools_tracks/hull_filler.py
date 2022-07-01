from starter2 import *
import convex_hull_tools as CHT
reload(CHT)
plt.close('all')
import three_loopers_six as TL
sim_list=['u601','u602','u603']
if 'ht' not in dir() :
    ht = {}
    for this_simname in sim_list:
        ht[this_simname] = CHT.hull_tool(TL.loops[this_simname])
        ht[this_simname].make_hulls() 
        ht[this_simname].make_overlaps()

if 'volumes' not in dir():
    volumes={}
    for sim in sim_list:
        volumes[sim]=CHT.filler( ht[sim],do_plots=False)#, core_list=[323]

if 0:
    for sim in volumes:
        fig,ax=plt.subplots(1,1)
        ax.scatter( volumes[sim]['volume_hull'], volumes[sim]['volume_fill'])
        ax.plot( volumes[sim]['volume_hull'], volumes[sim]['volume_hull'])
        fig.savefig('plots_to_sort/one_%s.png'%sim)
if 1:
    fig,ax = plt.subplots(1,1)

    for ns,sim in enumerate(sim_list):
        RRRR = nar(volumes[sim]['volume_part'])/nar(volumes[sim]['volume_hull'])
        RS = np.sort(RRRR)
        cuml = np.linspace(1,RRRR.size,RRRR.size)/RRRR.size
        ax.plot( RS,cuml,color=colors.color[sim])

    ax.legend(loc=0)
    axbonk(ax,xlabel=r'$V_{\rm{cell}}/V_{\rm{hull}}$',ylabel='N',yscale='log', xscale='linear')
    fig.savefig('plots_to_sort/cumlratios.png')

if 1:
    fig,ax = plt.subplots(1,2)
    ax0=ax[0];ax1=ax[1]

    for ns,sim in enumerate(sim_list):
        RRRR = nar(volumes[sim]['volume_part'])/nar(volumes[sim]['volume_hull'])
        bins = np.linspace( 0, 1,32)
        #hist,bins=np.histogram( RRRR , bins=bins)
        #bc = 0.5*(bins[1:]+bins[:-1])
        #ax0.plot(bc,hist,color=colors.color[sim], label=r'$sim%s$'%(ns+1))
        ax0.hist( RRRR, bins=bins,color=colors.color[sim], label=r'$sim%s$'%(ns+1), histtype='step')
        ax1.hist( nar(volumes[sim]['volume_fill'])/nar(volumes[sim]['volume_hull']),color=colors.color[sim], label=r'$sim%s$'%(ns+1), 
                histtype='step',bins=bins)

    ax0.legend(loc=0)
    axbonk(ax0,xlabel=r'$V_{\rm{cell,i}}/V_{\rm{hull,i}}$',ylabel='N',yscale='linear', xscale='linear')
    axbonk(ax1,xlabel=r'$V_{\rm{cell,all}}/V_{\rm{hull,i}}$',ylabel='N',yscale='linear', xscale='linear')
    ax1.set_ylim( ax0.get_ylim())
    fig.savefig('plots_to_sort/ratios.png')




