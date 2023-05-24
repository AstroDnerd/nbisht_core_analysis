from starter2 import *
import xtra_energy
import three_loopers_u500 as TL
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import pcolormesh_helper as pch
import colors
import movie_frames 
G = colors.G
plt.close('all')
def plotsa(accreter, tsing):


    fig, ax=plt.subplots(1,1)

    for core_id in accreter['rates']:
        ax.plot( accreter['time']/tsing.tsing_core[core_id], accreter['rates'][core_id])
    fig.savefig('plots_to_sort/doop.png')
    plt.close(fig)

    fig, ax=plt.subplots(1,1)
    #ax0=ax[0]; ax1=ax[1]
    ax0=ax

    for core_id in accreter['rates']:
        ax0.plot( accreter['time']/tsing.tsing_core[core_id], accreter['r10k'][core_id], c=[0.5]*4)
        #ax0.plot( accreter['time'], accreter['r10k'][core_id])
    ax0.set(yscale='log',ylabel=r'$R_{1000} [\rm{AU}]$', xlabel=r'$t/tsing$', xlim=[0,2])
    ax0.axvline(1,c='k',linewidth=0.1)
    m10k = nar(sorted(accreter['mass_10k']))
    cdf = np.arange(m10k.size)/m10k.size
    #ax1.plot( m10k, cdf)
    #ax1.hist(m10k)
    #ax1.set(xlabel='Mass[Msun]',ylabel='CDF')

    fig.savefig('plots_to_sort/mass_10k.pdf')



def accretion(things, tsing):



    rates={}
    r10k = {}
    mass_r10k=[]
    for core_id in things:
        #core_id=obj.core_id
        obj = things[core_id]
        #ax6.set(xlabel='R',xscale='log',ylabel='|V_r|',yscale='log')
        #ax6.set_title('Slope at tend: %0.2f'%pfit[0])
        #ax6.plot( r_bins[1:], 4*np.pi/3*r_bins[1:]**3,'k')
        dMdTf = obj.dMdTf+0
        MMM = obj.MMM+0
        ok = ~np.isnan(dMdTf)
        ext = extents(dMdTf[ok])
        #okt = ~np.isnan(MMM)
        #extt = extents(MMM[ok])
        #for i in range(3):
        #    ax6[i][0].set(ylim=ext.minmax)
        #    ax6[i][1].set(ylim=extt.minmax)
        #fig6.savefig('plots_to_sort/velocity_vs_radius_%s_c%04d.png'%(this_looper.sim_name, core_id))

        fig,axes=plt.subplots(3,1,figsize=(4,8))
        #fig,axes=plt.subplots(2,2,figsize=(12,12))
        #ax0=axes[0][0]; ax1=axes[0][1] ax2=axes[1][0]; ax3=axes[1][1]
        ax0=axes[0];ax1=axes[1];ax2=axes[2]
        fig.subplots_adjust(hspace=0)
        rcen = 0.5*(obj.r_bins[1:]+obj.r_bins[:-1])

        #
        # UNITS.  THERES BETTER PLACES FOR THIS.
        #
        #4.9 pc in AU
        rcen *= colors.length_units_au
        #5900 solar masses
        MMM *= colors.mass_units_msun

        times=obj.times.flatten()
        XXX,YYY = np.meshgrid( times,rcen)
        ok = ~np.isnan(MMM)
        radius_to_cut_off = 1000
        index = np.where( rcen < radius_to_cut_off)[0].max()
        mass_at_the_end = MMM[index,-1]
        mass_r10k.append(mass_at_the_end)
        fiducial=mass_at_the_end
        Max=MMM[ok].max()
        minner_chicken_dinner = fiducial**2/Max
        norm_mass = mpl.colors.LogNorm(vmin=minner_chicken_dinner,vmax=Max)

        #print(thalf.shape)
        #print(obj.times.shape)
        #print(MMM.shape)
        thalf = 0.5*(obj.times[1:]+obj.times[:-1])
        dt    =     (obj.times[1:]-obj.times[:-1])
        dt.shape=1,dt.size
        dMass = MMM[:,1:]-MMM[:,:-1]
        Mcen  = 0.5*(MMM[:,1:]+MMM[:,:-1])
        dMdt = dMass/dt/Mcen*colors.tff
        ok = ~np.isnan(dMdt)
        maxmax=np.abs(dMdt[ok]).max()
        maxmax=4
        norm2 = mpl.colors.Normalize(-maxmax,maxmax)

        


        ax=axes
        #norm = mpl.colors.LogNorm( dMdTf[dMdTf>0].min(), dMdTf.mean())
        ok = ~np.isnan(dMdTf)
        maxmax = np.abs(dMdTf[ok]).max()
        norm = mpl.colors.Normalize( -maxmax,maxmax)
        #pdb.set_trace()
        cmap=copy.copy(mpl.cm.get_cmap("seismic"))
        cmap.set_bad([0.9]*3)

        M3 = Mcen+0
        M3[ np.isnan(M3)]=1e9
        M10k = np.argmin(np.abs(M3-mass_at_the_end),axis=0)
        ind = np.arange(M3.shape[1])
        take = np.ravel_multi_index( nar([M10k,ind]), M3.shape)
        DDD = dMdt.flatten()[take]
        rates[core_id]=DDD
        ax2.plot(thalf.flatten(), DDD)

        #pdb.set_trace()
        print(M10k)
        ax1.plot(thalf.flatten(), rcen[M10k], c='k')
        r10k[core_id]=rcen[M10k]





        plot2=ax0.pcolormesh( XXX,YYY, MMM,  norm=norm_mass, shading='nearest', cmap=cmap)
        fig.colorbar(plot2, label=r'$M(<r)~~ [M_\odot]$',ax=ax0)

        AXX = ax1
        X2,Y2=np.meshgrid(thalf.flatten(),rcen)
        plot3=AXX.pcolormesh(X2 ,Y2, dMdt,  norm=norm2, shading='nearest', cmap=cmap)
        fig.colorbar(plot3,label=r'$\frac{d \ln M}{d t/tff}$', ax=AXX )

        if tsing:
            for aaa in axes.flatten():
                aaa.axvline(  tsing.tsing_core[core_id],c='k', linewidth=0.5)
                #ax1.axvline( tsing.tsing_core[core_id],c='k')
                aaa.axvline(  tsing.tend_core[core_id],c='k', linewidth=0.5)
                #ax1.axvline( tsing.tend_core[core_id],c='k')

        ax0.set(yscale='log', ylabel=r'$R~~[\rm{AU}]$', xticks=[])
        ax1.set(yscale='log',ylabel=r'$R~~[\rm{AU}]$', xticks=[], xlim=ax0.get_xlim())
        ax2.set(yscale='linear',ylabel=r'$d ln M /dt/tff$',xlabel='t/tff', xlim=ax0.get_xlim())#,ylabel='R[AU]')
        ax2.set_aspect(ax0.get_aspect())

        fig.tight_layout()
        fig.subplots_adjust(hspace=0)
        fig.savefig('plots_to_sort/accretion_%s_c%04d.pdf'%(obj.sim_name, core_id))

    output={'rates':rates, 'time':thalf.flatten()}
    output['radius_to_cut_off']=1000
    output['mass_10k']=mass_r10k
    output['r10k']=r10k
    return output
def plotmancer(obj, tsing):
        core_id=obj.core_id
        #ax6.set(xlabel='R',xscale='log',ylabel='|V_r|',yscale='log')
        #ax6.set_title('Slope at tend: %0.2f'%pfit[0])
        #ax6.plot( r_bins[1:], 4*np.pi/3*r_bins[1:]**3,'k')
        dMdTf = obj.dMdTf+0
        MMM = obj.MMM+0
        ok = ~np.isnan(dMdTf)
        ext = extents(dMdTf[ok])
        #okt = ~np.isnan(MMM)
        #extt = extents(MMM[ok])
        #for i in range(3):
        #    ax6[i][0].set(ylim=ext.minmax)
        #    ax6[i][1].set(ylim=extt.minmax)
        #fig6.savefig('plots_to_sort/velocity_vs_radius_%s_c%04d.png'%(this_looper.sim_name, core_id))



        fig,axes=plt.subplots(3,1,figsize=(4,8))
        #fig,axes=plt.subplots(2,2,figsize=(12,12))
        #ax0=axes[0][0]; ax1=axes[0][1] ax2=axes[1][0]; ax3=axes[1][1]
        ax0=axes[0];ax1=axes[1];ax2=axes[2]
        fig.subplots_adjust(hspace=0)
        rcen = 0.5*(obj.r_bins[1:]+obj.r_bins[:-1])

        #
        # UNITS.  THERES BETTER PLACES FOR THIS.
        #
        #4.9 pc in AU
        rcen *= colors.length_units_au
        #5900 solar masses
        MMM *= colors.mass_units_msun

        XXX,YYY = np.meshgrid( obj.times.flatten(),rcen)
        ax=axes
        #norm = mpl.colors.LogNorm( dMdTf[dMdTf>0].min(), dMdTf.mean())
        ok = ~np.isnan(dMdTf)
        maxmax = np.abs(dMdTf[ok]).max()
        norm = mpl.colors.Normalize( -maxmax,maxmax)
        #pdb.set_trace()
        cmap=copy.copy(mpl.cm.get_cmap("seismic"))
        cmap.set_bad([0.9]*3)

        ok = ~np.isnan(MMM)
        radius_to_cut_off = 1000
        index = np.where( rcen < radius_to_cut_off)[0].max()
        mass_at_the_end = MMM[index,-1]
        fiducial=mass_at_the_end
        Max=MMM[ok].max()
        minner_chicken_dinner = fiducial**2/Max
        norm_mass = mpl.colors.LogNorm(vmin=minner_chicken_dinner,vmax=Max)


        plot2=ax0.pcolormesh( XXX,YYY, MMM,  norm=norm_mass, shading='nearest', cmap=cmap)
        fig.colorbar(plot2, label=r'$M(<r)~~ [M_\odot]$',ax=ax0)
        if 0:
            #MASS FLUX FROM VELOCITY.
            #UNDERSTAND ME LATER.
            plot=ax1.pcolormesh( XXX,YYY, dMdTf,  norm=norm, shading='nearest', cmap=cmap)
            fig.colorbar(plot, label=r'$\frac{\dot{M}}{M/ t_{ff}}$',ax=ax1)
            ax1.set(yscale='log', ylabel = r'$R~~[\rm{AU}]$')

        AXX = ax1
        thalf = 0.5*(obj.times[1:]+obj.times[:-1])
        dt    =     (obj.times[1:]-obj.times[:-1])
        #print(thalf.shape)
        #print(obj.times.shape)
        #print(MMM.shape)
        dt.shape=1,dt.size
        dMass = MMM[:,1:]-MMM[:,:-1]
        Mcen  = 0.5*(MMM[:,1:]+MMM[:,:-1])
        dMdt = dMass/dt/Mcen*colors.tff
        ok = ~np.isnan(dMdt)
        maxmax=np.abs(dMdt[ok]).max()
        maxmax=4
        norm2 = mpl.colors.Normalize(-maxmax,maxmax)
        X2,Y2=np.meshgrid(thalf.flatten(),rcen)
        plot3=AXX.pcolormesh(X2 ,Y2, dMdt,  norm=norm2, shading='nearest', cmap=cmap)
        fig.colorbar(plot3,label=r'$\frac{d \ln M}{d t/tff}$', ax=AXX )

        Rcen=0.5*(dMdTf[:,1:]+dMdTf[:,:-1])
        fig6,ax666=plt.subplots(1,1)
        #ax666.scatter(Rcen, dMdt)
        ax666.hist(np.log(np.abs(Rcen.flatten())),histtype='step')
        ax666.hist(np.log(np.abs(dMdt.flatten())*0.10),histtype='step')
        print(dt)
        fig6.savefig('plots_to_sort/fork_c%04d.png'%core_id)


        newcmp = mpl.colors.LinearSegmentedColormap.from_list('custom blue', 
                                             [(0,         '#0000ff'),
                                              (norm(0.5), '#dddddd'),
                                              (norm(2.0), '#00ffff'),
                                              (1,         '#ff0000')], N=256)

        norm_grav = mpl.colors.LogNorm( vmin=0.01,vmax=100)
        ok = ~np.isnan(obj.RVg)
        AXX=ax2
        plot=AXX.pcolormesh( XXX,YYY, np.abs(obj.RVg),  norm=norm_grav, shading='nearest', cmap=cmap)
        fig.colorbar(plot,label='EG/EK', ax=AXX )
        ax0.set(yscale='log', ylabel=r'$R~~[\rm{AU}]$', xticks=[])
        ax2.set(yscale='log',ylabel=r'$R~~[\rm{AU}]$',xlabel='t/tff', xlim=ax0.get_xlim())#,ylabel='R[AU]')
        ax1.set(yscale='log',ylabel=r'$R~~[\rm{AU}]$', xticks=[], xlim=ax0.get_xlim())
        #pdb.set_trace()
        
        if 0:
            rrrr = msR.transpose()
            rrrr = rrrr[mask,:]

            ax.plot(times , rrrr, c=[0.5]*3, linewidth=0.1, alpha=0.5)
        if tsing:
            for aaa in axes.flatten():
                aaa.axvline(  tsing.tsing_core[core_id],c='k', linewidth=0.5)
                #ax1.axvline( tsing.tsing_core[core_id],c='k')
                aaa.axvline(  tsing.tend_core[core_id],c='k', linewidth=0.5)
                #ax1.axvline( tsing.tend_core[core_id],c='k')


        fig.tight_layout()
        fig.subplots_adjust(hspace=0)
        fig.savefig('plots_to_sort/mass_flux_color_%s_c%04d.pdf'%(obj.sim_name, core_id))
class Relax2():
    def __init__(self):
        pass
def thing_saver(things,sim):
    fptr = h5py.File('this_mass_flux_%s.h5'%sim,'w')
    try:
        for core_id in things:
            obj = things[core_id]
            core_group=fptr.create_group('c%04d'%core_id)
            core_group['times']=obj.times
            core_group['core_id']=obj.core_id
            core_group['r_bins']=obj.r_bins
            core_group['dMdTf']=obj.dMdTf
            core_group['MMM']=obj.MMM
            core_group['RVg']=obj.RVg
            core_group['divv']=obj.divv
            core_group['meanvr']=obj.meanvr
            core_group['sim_name']=obj.sim_name
    except:
        raise
    finally:
        fptr.close()


def thing_reader(fname):
    fptr=h5py.File(fname,'r')
    bucket={}
    try:
        for group in fptr:
            this_one = Relax2()
            for quan in fptr[group]:
                QQQ = fptr[group][quan]
                if quan == 'sim_name':
                    QQQ = QQQ.asstr()
                this_one.__dict__[quan]=QQQ[()]
            bucket[this_one.core_id] = this_one

            
    except:
        raise
    finally:
        fptr.close()
    return bucket
class Relaxor():
    def __init__(self,this_looper):
        self.this_looper=this_looper
    def run(self,core_list=None, do_plots=True, r_inflection=None, frame_list=None, tsing=None):
        this_looper=self.this_looper
        self.sim_name = this_looper.sim_name

        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        if frame_list is None:
            thtr=this_looper.tr
            mask = movie_frames.quantized_mask(this_looper).flatten()
            ok = np.zeros_like(mask)
            ok[::10] = mask[::10]
            #mask=ok ;print('kludge mask')
            times=thtr.times[mask]+0 #the zero makes a copy
            times.shape=times.size,1
            times=times/colors.tff
            self.times=times
            frame_list=thtr.frames[mask]

        rm = rainbow_map(len(frame_list))
        for core_id in core_list:
            self.core_id=core_id
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            
            r_bins = np.geomspace( 2e-4, 32/128, 32)
            self.r_bins=r_bins
            r_cen = 0.5*(r_bins[1:]+r_bins[:-1])
            self.dMdTf = np.zeros( [len(r_bins)-1, len(frame_list)])
            self.dMdT3d = np.zeros( [len(r_bins)-1, len(frame_list)])
            self.MMM = np.zeros( [len(r_bins)-1, len(frame_list)])
            self.RVg = np.zeros( [len(r_bins)-1, len(frame_list)])
            self.divv=np.zeros_like(self.MMM)
            self.meanvr=np.zeros_like(self.MMM)
            frame_rmap=rainbow_map(len(frame_list))
            fig6,ax6=plt.subplots(3,2)
            ax_index=0
            got_tsing=False
            got_tend=False
            for nframe,frame in enumerate(frame_list):
                print("Flux on %s c%04d n%04d"%(this_looper.sim_name, core_id, frame))
                ds = this_looper.load(frame)
                xtra_energy.add_energies(ds)
                xtra_energy.add_gravity(ds)
                nf = np.where( this_looper.tr.frames == frame)[0][0]
                time = times[nframe,0]

                c = nar([ms.mean_xc[nf], ms.mean_yc[nf],ms.mean_zc[nf]])
                msR = ms.rc
                msR[ msR<1/2048]=1/2048
                
                MaxRadius=msR[:,nf].max()
                Radius = max([8.0/128, MaxRadius])
                rsph = ds.arr(Radius,'code_length')
                sp = ds.sphere(c,rsph)

                import other_scrubber
                reload(other_scrubber)


                R1 = sp['radius']
                order = np.argsort(R1)
                vel = []
                for axis in 'xyz':
                    vel.append( sp['velocity_%s'%axis][order][:10].mean())
                scrub = other_scrubber.scrubber(sp, reference_velocity = vel)
                scrub.compute_ge()
                scrub.compute_ke_rel()


                #Get data arrays

                dv = scrub.cell_volume
                RR = scrub.r
                DD = scrub.density
                divv = sp['velocity_divergence']
                #dv = np.abs(sp[YT_cell_volume])
                #RR =sp[YT_radius]
                #DD = sp[YT_density]

                ORDER = np.argsort( RR)
                RR_cuml = RR[ORDER]
                M_local = DD[ORDER]*dv[ORDER]
                M_cuml = np.cumsum( DD[ORDER]*dv[ORDER])
                V_cuml = np.cumsum( dv[ORDER])
                V_local = dv[ORDER]

                EG = scrub.ge
                EK = scrub.ke_rel

                EG_cuml = np.cumsum( EG[ORDER]*V_local)/V_cuml
                EK_cuml = np.cumsum( EK[ORDER]*V_local)/V_cuml

                vr = scrub.vr_rel
                flux = -( DD[ORDER]*vr[ORDER]*dv[ORDER])
                #flux = -np.cumsum( DD[ORDER]*vr[ORDER]*dv[ORDER])


                digitized = np.digitize( RR_cuml, r_bins)
                mean_flux  =nar([ flux[ digitized == i].sum() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])
                mass_quant  =nar([ M_cuml[ digitized == i].mean() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])
                EGmeans =nar([ EG_cuml[ digitized == i].mean() if (digitized==i).any() else -4000 for i in range(1,len(r_bins))])
                EKmeans =nar([ EK_cuml[ digitized == i].mean() if (digitized==i).any() else -4000 for i in range(1,len(r_bins))])



                divv_cuml = np.cumsum(divv[ORDER]*V_local)
                total_vr  =nar([ vr[ digitized == i].sum() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])
                total_divv  =nar([ divv_cuml[ digitized == i].mean() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])

                mass_flux_3d = np.cumsum(-(DD[ORDER]*divv[ORDER]*dv[ORDER]))
                F_3d  =nar([ mass_flux_3d[ digitized == i].mean() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])
                ok = ~np.isnan(F_3d)
                self.dMdT3d[ok,nframe] = F_3d[ok]/mass_quant[ok]*colors.tff

                ok = ~np.isnan(total_vr)
                self.meanvr[ok,nframe]=total_vr[ok]
                self.meanvr[~ok,nframe]=np.nan
                ok = ~np.isnan(total_divv)
                self.divv[ok,nframe]=total_divv[ok]
                self.divv[~ok,nframe]=np.nan



                ok = ~np.isnan(mean_flux)
                self.dMdTf[ok,nframe]=mean_flux[ok]/mass_quant[ok]*colors.tff
                self.dMdTf[~ok,nframe]=np.nan

                ok = ~np.isnan(mass_quant)
                self.MMM[ok,nframe]=mass_quant[ok]
                self.MMM[~ok,nframe]=np.nan

                ok = EKmeans>0
                self.RVg[ok,nframe]=(EGmeans/EKmeans)[ok]
                self.RVg[~ok,nframe]=np.nan

                #r_bins_c = 0.5*(r_bins[1:]+r_bins[:-1])
                #ax6.plot(r_bins_c, vr_mean, c = frame_rmap(nframe))
                if 0:
                    linewidth=0.3
                    c=frame_rmap(nframe)
                    if tsing:
                        if not got_tsing:
                            if time > tsing.tsing_core[core_id]:
                                c = 'g'
                                linewidth=1
                                got_tsing=True
                                ax_index=1
                        if not got_tend:
                            if time > tsing.tend_core[core_id]:
                                c = 'r'
                                linewidth=1
                                got_tend=True
                                pfit = np.polyfit( np.log10(RR_cuml[1:]), np.log10(M_cuml[1:]/RR_cuml[1:]**3),1)
                                print(pfit)
                                ax_index=2

                    from scipy.ndimage import gaussian_filter

                    ##ax6.plot(RR_cuml, gaussian_filter(M_cuml,1), c=c,linewidth=linewidth)
                    ##ax6.plot(RR_cuml[1:], M_cuml[1:]/RR_cuml[1:]**3, c=c,linewidth=linewidth)
                    ax6[ax_index][0].plot( r_cen[ok], vr_mean[ok], c=c, linewidth=linewidth)
                    ax6[ax_index][1].plot( r_cen[ok], vt_mean[ok], c=c, linewidth=linewidth)
                    #ax6[0].hist( RR_cuml.v, histtype='step', color=c)
            plotmancer(self,tsing)



if 0:
    import three_loopers_six as TL
    sim_list=['u601','u602','u603']
    sim_list=['u602']
import three_loopers_u500 as TL
sim_list=['u501','u502','u503']
sim_list=['u502']

import tsing
reload(tsing)
if 'tsing_tool' not in dir():
    tsing_tool={}
    for ns,sim in enumerate(sim_list):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()

#import anne
#reload(anne)
##anne.make_inflection()
if 'things' not in dir():
    things={}

if 1:
    #save and replot.
    if 0:
        thing_saver(things,sim)
    if 1:
        thing_new = thing_reader('this_mass_flux_u502.h5')
    if 0:
        sim='u502'
        accreter=accretion(thing_new, tsing_tool[sim])
    if 1:
        plotsa(accreter,tsing_tool[sim])

    if 1:
        #paper plot.  Read off disk.
        plotmancer( thing_new[114], tsing_tool[sim])

if 0:
    #new plots.
    for sim in sim_list:
        #all_cores=np.unique( TL.loops[sim].tr.core_ids)
        #core_list=list(all_cores)
        #core_list=None#[323]
        #core_list=[323]
        #core_list=[25]
        #core_list=[74]
        #core_list=[114]
        #core_list=[195]
        #core_list=[361]
        #core_list=[8]
        #core_list=[381]
        #core_list=[323]
        #core_list=None
        #core_list = [9]
        core_list = TL.loops[sim].core_by_mode['Alone']
        #core_list=core_list[4:7]

        core_list = [114]
        for core_id in core_list:
            if core_id in things:
                continue
            thing = Relaxor( TL.loops[sim])
            thing.run(core_list=[core_id], do_plots=True, frame_list=None, tsing=tsing_tool[sim])#, r_inflection=anne.inflection[sim])
            things[core_id]=thing

        plotmancer(thing, tsing_tool[sim])

    if 1:
        thing_saver(things,sim)
    if 0:
        if 'accreter' not in dir() or True:
            accreter=accretion(things, tsing_tool[sim])
    if 0:
        plotsa(accreter,tsing_tool[sim])

    if 0:
        for core_id in things:
            #plotmancer(things[core_id], tsing_tool[sim])
            calculus(things[core_id],tsing_tool[sim])

if 0:
    plotmancer(thing,tsing_tool[sim])
if 0:
    accretion(thing,tsing_tool[sim])
