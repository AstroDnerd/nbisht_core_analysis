from starter2 import *
from collections import defaultdict
import heat_map
from scipy.optimize import curve_fit
import colors
reload(heat_map)
def b_from_rho(lnrho,slope,offset):
    return slope*lnrho+offset
class brho():
    def __init__(self, loop):
        self.this_looper=loop
        self.cores_used=[]
        self.mean_b=defaultdict(list)
        self.mean_rho=defaultdict(list)
        self.fit1=defaultdict(list)
        self.alpha_FTA=defaultdict(list)
        self.alpha_ATF=[]
        
    def run(self,core_list=None,plot=False):
        dx=1./2048
        nx = 1./dx
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        self.times=tsorted

        self.b_array=np.zeros([len(core_list), len(tsorted)])
        self.d_array=np.zeros([len(core_list), len(tsorted)])
        self.a_array=np.zeros([len(core_list), len(tsorted)])
        for nc,core_id in enumerate(core_list):
            print('do %s %d'%(self.this_looper.sim_name,core_id))
            #continue    
            #ms = trackage.mini_scrubber(thtr,core_id)
            #ms.particle_pos(core_id)
            #self.ms = ms
            #if ms.nparticles < 10:
            #    continue

            self.cores_used.append(core_id)

            B = thtr.c([core_id],'magnetic_field_strength')
            density = thtr.c([core_id],'density')
            cell_volume = thtr.c([core_id],'cell_volume')
            Vtotal=cell_volume.sum(axis=0)
            self.mean_b[core_id] =  (B*cell_volume).sum(axis=0)/Vtotal
            self.b_array[nc,:]=self.mean_b[core_id]
            self.mean_rho[core_id]= (density*cell_volume).sum(axis=0)/Vtotal
            self.d_array[nc,:]=self.mean_rho[core_id]
            rm=rainbow_map(tsorted.size)
            if plot:
                fig,ax=plt.subplots(1,1)

            for nframe in range(B.shape[1]-1):
                BBB=np.log10(B[:,nframe])
                DDD=np.log10(density[:,nframe])
                popt,pcov = curve_fit( b_from_rho, DDD, BBB, p0=[1,1])
                self.alpha_FTA[core_id].append(popt[0])

                if plot:

                    ax.scatter(DDD , BBB, c=[rm(nframe)]*density.shape[0])
                    #ax.plot(DDD, DDD*pfit[1]+pfit[0])
                    ax.plot( DDD, b_from_rho(DDD,*popt),c='k')
            if np.isnan(self.alpha_FTA[core_id]).any():
                pdb.set_trace()
            self.a_array[nc,:-1]=self.alpha_FTA[core_id]
            if plot:
                fig.savefig('plots_to_sort/b_rho_%s_c%04d.png'%(self.this_looper.sim_name, core_id))



if 0:
    import three_loopers_six as TL6
    thing = flow(TL6.loops['u603'])
    core_list=[76]
    thing.run(core_list=core_list)

if 'brthings' not in dir():
    import three_loopers_u500 as TL5
    simlist = ['u501','u502','u503']
    brthings={}
    for sim in simlist:
        brthings[sim]= brho(TL5.loops[sim])
        brthings[sim].run()


if 0:
    import three_loopers_six as TL6
    brthing = brho(TL6.loops['u603'])
    brthing.run()

if 0:
    #make ATF
    fig,ax=plt.subplots(1,1)
    d_array=brthing.d_array
    b_array=brthing.b_array
    a_array=brthing.a_array

    alpha_ATF=[]
    for nf in np.arange(0,120,10):
        ddd = np.log10(d_array[:,nf])
        bbb = np.log10(b_array[:,nf])
        aaa = a_array[:,nf]
        ax.scatter( ddd,bbb)
        drho=0.02*np.mean(ddd)
        d1 = ddd-drho
        d2 = ddd+drho
        b1=bbb-aaa*drho
        b2=bbb+aaa*drho
        dlines=np.stack([d1,d2])
        blines=np.stack([b1,b2])
        plt.plot(dlines,blines)
        outname='plots_to_sort/b_rho_WTF_%s_i%04d.png'%(this_looper.sim_name,nf)
        fig.savefig(outname)
        print(outname)

if 0:
    #mean B
    import heat_map
    reload(heat_map)
    fig,ax=plt.subplots(1,1)
    this_looper=brthing.this_looper
    B = this_looper.tr.track_dict['magnetic_field_strength']
    minmin = B.min()
    maxmax = B.max()
    bins = np.geomspace(minmin,maxmax,64)
    heat_map.plot_heat(times = this_looper.tr.times, cores_used = brthing.cores_used, quan_dict=brthing.mean_b ,ax=ax, cut_last=True,bins=bins)
    axbonk(ax,yscale='log',xlabel='T', ylim=[minmin,maxmax])
    fig.savefig('plots_to_sort/mean_b_per_core_%s.png'%brthing.this_looper.sim_name)

if 0:
    #make ATF
    atf={}
    for sim in brthings:
        brthing  = brthings[sim]
        fig,ax=plt.subplots(1,1)
        d_array=brthing.d_array
        b_array=brthing.b_array

        atf[sim] = []
        for nt in range(d_array.shape[1]):
            DDD = np.log10(d_array[:,nt])
            BBB = np.log10(b_array[:,nt])
            popt,pcov = curve_fit( b_from_rho, DDD, BBB, p0=[1,1])

            atf[sim].append(popt[0])
            ax.scatter( DDD,BBB)
            ax.plot( DDD, b_from_rho(DDD,*popt))
        fig.savefig('plots_to_sort/b_rho_AFT_%s.png'%sim)

if 1:
    #plot FTA and overplot ATF.
    import heat_map
    reload(heat_map)
    fig,axes=plt.subplots(1,3, figsize=(12,8), sharey=True)
    axlist=axes.flatten()
    for ns,sim in enumerate(brthings):
        ax=axlist[ns]
        brthing=brthings[sim]
        this_looper=brthing.this_looper
        heat_map.plot_heat(times = this_looper.tr.times[:-1], cores_used = brthing.cores_used , quan_dict=brthing.alpha_FTA,ax=ax, cut_last=True)
        ax.plot((this_looper.tr.times/colors.tff)[:-1], atf[sim][:-1],c='k')

    fig.savefig('plots_to_sort/alpha_per_core_%s.png'%brthing.this_looper.sim_name)

