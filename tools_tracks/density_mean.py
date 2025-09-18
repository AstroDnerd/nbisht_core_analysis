
from starter2 import *
import data_locations as dl
from collections import defaultdict

import davetools
reload(davetools)
import mass_tools
reload(mass_tools)

if 1:
    import three_loopers_mountain_top as TLM
    reload(TLM)
if 'rt_dict' not in dir() :
    rt_dict={}
    simlist = ['u301']#,'u302','u303']
    for this_simname in simlist:

        rt_dict[this_simname] = mass_tools.mass_tool(TLM.loops[this_simname])
        rt_dict[this_simname].run()#core_list=TLM.loops[this_simname].core_list[:10])
if 1:
    for this_simname in simlist:
        fig,ax=plt.subplots(2,2,figsize=(12,12))

        mt = rt_dict[this_simname]

        volume = np.zeros([ len(mt.cores_used), len(mt.times)])
        density = np.zeros_like(volume)
        for nt,core_id in enumerate( mt.cores_used):
            volume[nt,:]=mt.volume[core_id]
            density[nt,:]=mt.mean_rho[core_id]

        norm = mpl.colors.Normalize()
        n0=0
        norm.autoscale( np.log10(density[:,n0]))
        cmap = mpl.cm.jet
        color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)

        mean_density = 10**( np.log10( density).mean(axis=0))
        density_cmap = rainbow_map( np.log10(density)[:,0].max())
        mass_ratio = 1
        big_core = 0
        core_list = []
        for nc,core_id in enumerate(mt.cores_used):
            ax[0][0].plot(mt.times, mt.volume[core_id],c=[0.5]*4)
            #c=density_cmap( np.log10(mt.mean_rho[core_id][0]))
            c = color_map.to_rgba(np.log10(density[nc,n0]))
            Y2=mt.mean_rho[core_id]#/nar(mt.mean_rho[core_id])[0:2].mean()
            ax[0][1].plot(mt.times, Y2,c=c)#[0.5]*4)

            M = mt.unique_mass[core_id] #/nar(mt.unique_mass[core_id])[:6].mean()
            Mbar = mt.unique_mass[core_id]/nar(mt.unique_mass[core_id])[:6].mean()
            this_ratio = Mbar.max()/Mbar.min()
            if Mbar.min() < 0.1:
                core_list.append(nc)

            ax[1][0].plot(mt.times,M)
        mass_tools.plot_mass_tracks(mt,ax[1][1],core_list=core_list)


            #ax[1].plot([mt.times[0], mt.times[-1]], [Y2[0],Y2[-1]],c=c)#[0.5]*4)
            #ax[1].plot(mt.times, mt.mean_rho_w[core_id],c=[0.5]*4)

        axbonk(ax[1][0],xlabel='t',ylabel='M',yscale='log')
        axbonk(ax[0][1],xlabel='t',ylabel='rho',yscale='log')
        axbonk(ax[0][0],xlabel='t',ylabel='V',yscale='log')
        fig.savefig('plots_to_sort/%s_meansies.pdf'%this_simname)
        plt.close(fig)
