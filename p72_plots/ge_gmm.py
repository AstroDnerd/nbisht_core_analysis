
from starter2 import *
import xtra_energy
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter

from collections import defaultdict

#import three_loopers_u500 as TL
import three_loopers_six as TL
plt.close('all')


class ge_tool():
    def __init__(self, this_looper, inflection):
        self.this_looper=this_looper
        self.output = defaultdict(list)
        self.rinflection = inflection
    def run(self,core_list=None):
        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        this_looper=self.this_looper
        frame = this_looper.target_frame
        ds = this_looper.load(frame)
        G = ds['GravitationalConstant']/(4*np.pi)
        xtra_energy.add_energies(ds)
        for core_id in core_list:
            print('GE tool %s %d'%(this_looper.sim_name,core_id))

            ms = trackage.mini_scrubber(this_looper.tr,core_id)

            c = nar([ms.mean_x[-1], ms.mean_y[-1],ms.mean_z[-1]])
            
            R_SPHERE = 8/128
            rsph = ds.arr(R_SPHERE,'code_length')
            sp = ds.sphere(c,rsph)

            GE = np.abs(sp['grav_energy'])
            GE2 = np.abs(sp[YT_grav_energy_2])
            dv = sp['cell_volume']
            RR = sp['radius']
            DD = sp['density']

            R_KEEP = self.rinflection[core_id]

            #get the zones in the sphere that are
            #within R_KEEP
            ok_fit = np.logical_and(RR  < R_KEEP, GE>0)

            rok=RR[ok_fit].v
            #rok[rok<1/2048] = 1/2048

            #k1
            #ok_fit = np.logical_and(ok_fit, RR  > r_cen[keepers].min())

            ge_total = (GE[ok_fit]*dv[ok_fit]).sum()
            ge_total_2 = (GE2[ok_fit]*dv[ok_fit]).sum()
            mtotal = (sp['cell_mass'][ok_fit]).sum()
            gmm      = G*mtotal**2/R_KEEP
            self.output['ge_total'].append(ge_total)
            self.output['ge_total_2'].append(ge_total_2)
            self.output['gmm'].append(gmm)

            if 0:
                #Better ge fit
                def plain_powerlaw( x, q, norm):
                    return (2*q+2)*x+norm#np.log10(G*mtotal/R_KEEP)
                popt, pcov=curve_fit(plain_powerlaw, np.log10(rok/R_KEEP), np.log10(GE[ok_fit]))
                A = 10**popt[1]
                alpha = popt[0]
                self.output['alpha_ge'].append(alpha)
                self.output['A'].append(A)

sim_list=['u601']
if 1:
    import three_loopers_six as TL
    import r_inflection
    reload(r_inflection)
    if 'inflection' not in dir():
        inflection = {}
    for sim in sim_list:
        if sim in inflection:
            continue
        inflection[sim]=r_inflection.R_INFLECTION( TL.loops[sim])
        inflection[sim].run()

if 'stuff' not in dir():
    stuff={}
    for sim in sim_list:
        #stuff={}
        all_cores=np.unique( TL.loops[sim].tr.core_ids)
        #core_list=list(all_cores)
        core_list=all_cores
        #core_list=[323]
        #stuff[sim]=plot_phi( TL.loops[sim],core_list=core_list, do_plots=False)
        stuff[sim]=ge_tool(TL.loops[sim], inflection[sim].rinflection)
        stuff[sim].run(core_list=core_list)

if 0:
    fig,ax=plt.subplots(1,2)
    axlist=ax.flatten()

    for ns,sim in enumerate(sim_list):
        ss = stuff[sim]
        axlist[ns].scatter( ss.output['ge_total'], ss.output['gmm'])
        axlist[ns].plot( ss.output['ge_total'],ss.output['ge_total'])
        axlist[1].scatter( ss.output['ge_total'], nar(ss.output['gmm'])/nar(ss.output['ge_total']))
        axbonk(axlist[1],xscale='log')
        axbonk(axlist[ns],xlabel=r'$(\nabla \phi)^2$', ylabel='GMM/R',xscale='log',yscale='log')
    fig.savefig('plots_to_sort/ge_gmm.pdf')

if 1:
    fig,ax=plt.subplots(1,1)
    ax.scatter( inflection[sim].rinflection_list, nar(ss.output['ge_total'])/nar(ss.output['gmm']))
    axbonk(ax,xscale='log',yscale='log',xlabel='Rinf',ylabel='GMM/GE')
    fig.savefig('plots_to_sort/ge_radius.pdf')

if 0:
    #good.  Things changed but not a ton.
    fig,ax=plt.subplots(1,1)
    ax.hist( nar(ss.output['gmm'])/nar(ss.output['ge_total']), histtype='step', color='g')
    ax.hist( nar(gmm)/nar(ge), histtype='step',color='r')
    fig.savefig('plots_to_sort/gegegege.pdf')
    

if 0:
    #very good.  New gravity is good, not plagued by grid boundaries
    fig,ax=plt.subplots(1,2)
    axlist=ax.flatten()

    for ns,sim in enumerate(sim_list):
        ss = stuff[sim]
        axlist[ns].scatter( ss.output['ge_total'], nar(ss.output['ge_total'])/nar(ss.output['ge_total_2']))
        axbonk(axlist[ns],xlabel=r'$(\nabla \phi)^2$', ylabel='GMM/R',xscale='log',yscale='log')
    fig.savefig('plots_to_sort/ge_ge2.pdf')

