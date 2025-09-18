from starter2 import *
from collections import defaultdict
import heat_map
from scipy.optimize import curve_fit
import colors
reload(heat_map)
import pcolormesh_helper as pch
reload(pch)
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
        if plot:
            ds = this_looper.load(frame)
            ad = ds.all_data()
            b = ad[YT_magnetic_field_strength].v
            d = ad[YT_density].v
            dv= ad[YT_cell_volume].v
            dbins = np.geomspace( d[d>0].min(), d.max(),64)
            bbins = np.geomspace( b[b>0].min(), b.max(),64)
            hist, xb, yb = np.histogram2d( d, b, bins=[dbins,bbins],weights=dv)

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

                pch.helper(hist,xb,yb,ax=ax)

            for nframe in [0]:#range(B.shape[1]-1):
                BBB=np.log10(B[:,nframe])
                DDD=np.log10(density[:,nframe])
                popt,pcov = curve_fit( b_from_rho, DDD, BBB, p0=[1,1])
                self.alpha_FTA[core_id].append(popt[0])

                if plot:
                    ax.scatter(10**DDD , 10**BBB, c=[rm(nframe)]*density.shape[0])
                    #ax.scatter(DDD , BBB, c=[rm(nframe)]*density.shape[0])
                    #ax.plot(DDD, DDD*pfit[1]+pfit[0])
                    ax.plot( 10**DDD, 10**b_from_rho(DDD,*popt),c='k')

            if np.isnan(self.alpha_FTA[core_id]).any():
                pdb.set_trace()
            self.a_array[nc,:-1]=self.alpha_FTA[core_id]
            if plot:
                axbonk(ax,xscale='log',yscale='log',xlabel='density',ylabel='magnetic_field_strength')
                fig.savefig('plots_to_sort/b_rho_%s_c%04d.png'%(self.this_looper.sim_name, core_id))



if 0:
    import three_loopers_six as TL6
    thing = flow(TL6.loops['u603'])
    core_list=[76]
    thing.run(core_list=core_list)

import three_loopers_u500 as TL5
simlist = ['u501','u502','u503']
simlist = ['u502']
if 'brthings' not in dir() or True:
    brthings={}
if 1:
    for sim in simlist:
        if sim in brthings:
            continue
        core_list= TL5.loops[sim].core_by_mode['Alone']

        temp_thing= brho(TL5.loops[sim])
        #brthings[sim].run(core_list=[32],plot=True)
        #brthings[sim].run(plot=True)
        temp_thing.run(plot=True, core_list=core_list)
        brthings[sim]=temp_thing

if 0:
    nf=0
    frame=0
    farts=[]
    b_ext = extents()
    d_ext = extents()
    for nsim,sim_name in enumerate(simlist):
        brthing = brthings[sim_name]
        this_looper = brthings[sim_name].this_looper
        ds = this_looper.load(frame)
        ad = ds.all_data()
        b = ad[YT_magnetic_field_strength].v
        b_ext(b)
        d = ad[YT_density].v 
        d_ext(d)
    for nsim,sim_name in enumerate(simlist):
        fig,axes=plt.subplots(1,2)
        ax=axes[0]; ax1=axes[1]
        brthing = brthings[sim_name]
        this_looper = brthings[sim_name].this_looper
        print(sim_name)
        ds = this_looper.load(frame)
        print(ds.fullpath)
        ad = ds.all_data()
        b = ad[YT_magnetic_field_strength].v
        d = ad[YT_density].v 
        dv= ad[YT_cell_volume].v
        #dbins = np.geomspace( d[d>0].min(), d.max(),64)
        #bbins = np.geomspace( b[b>0].min(), b.max(),65)
        dbins = np.geomspace( d_ext.minmax[0], d_ext.minmax[1],64)
        bbins = np.geomspace( b_ext.minmax[0], b_ext.minmax[1],65)


        hist, xb, yb = np.histogram2d( d, b, bins=[dbins,bbins],weights=dv)
        #pch.helper(hist,xb,yb,ax=ax)
        
        d = this_looper.tr.track_dict['density'][:,nf]
        b = this_looper.tr.track_dict['magnetic_field_strength'][:,nf]
        dv = this_looper.tr.track_dict['cell_volume'][:,nf]
        hist, xb, yb = np.histogram2d( d, b, bins=[dbins,bbins],weights=dv)
        pch.helper(hist,xb,yb,ax=ax1)
        axbonk(ax,xscale='log',yscale='log',xlabel='density',ylabel='magnetic', xlim=d_ext.minmax,ylim=b_ext.minmax)
        axbonk(ax1,xscale='log',yscale='log',xlabel='density',ylabel='magnetic', xlim=d_ext.minmax,ylim=b_ext.minmax)
        #ax1.set_xlim( ax.get_xlim())
        #ax1.set_ylim( ax.get_ylim())

        fig.savefig('plots_to_sort/b_rho_ab_%s.png'%sim)


if 0:
    import three_loopers_six as TL6
    brthing = brho(TL6.loops['u603'])
    brthing.run()

if 0:
    #make ATF
    for sim_name in brthings:
        brthing = brthings[sim_name]
        this_looper = brthing.this_looper
        fig,ax=plt.subplots(1,1)
        d_array=brthing.d_array
        b_array=brthing.b_array
        a_array=brthing.a_array

        alpha_ATF=[]
        #frame_list = this_looper.tr.frames[::10]
        frame_list = this_looper.tr.frames[:1]
        
        for nf,frame in enumerate(frame_list):

            if 1:
                ds = this_looper.load(frame)
                ad = ds.all_data()
                b = ad[YT_magnetic_field_strength].v
                d = ad[YT_density].v
                dv= ad[YT_cell_volume].v
                dbins = np.geomspace( d[d>0].min(), d.max(),64)
                bbins = np.geomspace( b[b>0].min(), b.max(),64)
                hist, xb, yb = np.histogram2d( d, b, bins=[dbins,bbins],weights=dv)
                pch.helper(hist,xb,yb,ax=ax)
            #ddd = np.log10(d_array[:,nf])
            #bbb = np.log10(b_array[:,nf])
            ddd = d_array[:,nf]
            bbb = b_array[:,nf]
            aaa = a_array[:,nf]
            ax.scatter( ddd,bbb)
            drho=0.02*np.mean(ddd)
            d1 = ddd/drho
            d2 = ddd*drho
            b1=bbb/drho**aaa
            b2=bbb*drho**aaa
            dlines=np.stack([d1,d2])
            blines=np.stack([b1,b2])
            #plt.plot(dlines,blines)
            outname='plots_to_sort/b_rho_WTF_%s_i%04d.png'%(this_looper.sim_name,nf)
            axbonk(ax,xscale='log',yscale='log',xlabel='density',ylabel='magnetic field strength')
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

if 0:
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

