
from starter2 import *
import xtra_energy
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import pcolormesh_helper as pch
import colors
import movie_frames 
G = colors.G
plt.close('all')



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
        if len(core_list) != 1:
            print("Error: this tool only does one core at a time.")
            pdb.set_trace()
        for core_id in core_list:
            self.core_id=core_id
            ms = trackage.mini_scrubber(this_looper.tr,core_id,do_velocity=False)
            
            r_bins = np.geomspace( 2e-4, 32/128, 32)
            self.r_bins=r_bins
            r_cen = 0.5*(r_bins[1:]+r_bins[:-1])
            self.dMdTf = np.zeros( [len(r_bins)-1, len(frame_list)])
            self.dMdT3d = np.zeros( [len(r_bins)-1, len(frame_list)])
            self.MMM = np.zeros( [len(r_bins)-1, len(frame_list)])
            self.RVg = np.zeros( [len(r_bins)-1, len(frame_list)])
            self.Mcont={}
            self.Rcont={}
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

                self.Mcont[nframe]=M_cuml
                self.Rcont[nframe]=RR_cuml

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
            #plotmancer(self,tsing)


def plotmancer(obj, tsing=None, tall=True):
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

        if tall:
            nup = 3
            nacross = 1
            size=(4,8)
        else:
            nup=1
            nacross = 3
            size=(12,4)

        fig,axes=plt.subplots(nup,nacross,figsize=size)
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

        TIMES = obj.times
        XXX,YYY = np.meshgrid( obj.times.flatten(),rcen)
        ax=axes
        #norm = mpl.colors.LogNorm( dMdTf[dMdTf>0].min(), dMdTf.mean())
        ok = ~np.isnan(dMdTf)
        maxmax = np.abs(dMdTf[ok]).max()
        norm = mpl.colors.Normalize( -maxmax,maxmax)
        #pdb.set_trace()
        cmap=copy.copy(mpl.cm.get_cmap("seismic"))
        cmap.set_bad([0.9]*3)

        # in progress to detect early sim nan values...
        ok = ~np.isnan(MMM)
        radius_to_cut_off = 1000
        index = np.where( rcen < radius_to_cut_off)[0].max()
        min_mass_index = np.where(~np.isnan(MMM[index:,-1]))[0][0]+index
        index = max([min_mass_index,index])
        mass_at_the_end = MMM[index,-1]
        fiducial=mass_at_the_end
        Max=MMM[ok].max()
        minner_chicken_dinner = fiducial**2/Max
        #pdb.set_trace()
        norm_mass = mpl.colors.LogNorm(vmin=minner_chicken_dinner,vmax=Max)


        plot2=ax0.pcolormesh( XXX,YYY, MMM,  norm=norm_mass, shading='nearest', cmap=cmap)
        fig.colorbar(plot2, label=r'$M(<r)~~ [M_\odot]$',ax=ax0)
        if 0:
            #MASS FLUX FROM VELOCITY.
            #UNDERSTAND ME LATER - mass flux.
            plot=ax1.pcolormesh( XXX,YYY, dMdTf,  norm=norm, shading='nearest', cmap=cmap)
            fig.colorbar(plot, label=r'$\frac{\dot{M}}{M/ t_{ff}}$',ax=ax1)
            ax1.set(yscale='log', ylabel = r'$R~~[\rm{AU}]$')
        if 1:
            #dMdt
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
            maxmax=3.8
            norm2 = mpl.colors.Normalize(-maxmax,maxmax)
            X2,Y2=np.meshgrid(thalf.flatten(),rcen)
            plot3=AXX.pcolormesh(X2 ,Y2, dMdt,  norm=norm2, shading='nearest', cmap=cmap)
            fig.colorbar(plot3,label=r'$\frac{d \ln M}{d t/t_{\rm{ff}}}$', ax=AXX )
        if 0:
            #why is \rho v_r not Mdot?
            #understand me later. (same as above.)
            Rcen=0.5*(dMdTf[:,1:]+dMdTf[:,:-1])
            fig6,ax666=plt.subplots(1,1)
            #ax666.scatter(Rcen, dMdt)
            ax666.hist(np.log(np.abs(Rcen.flatten())),histtype='step')
            ax666.hist(np.log(np.abs(dMdt.flatten())*0.10),histtype='step')
            print(dt)
            fig6.savefig('plots_to_sort/fork_c%04d.png'%core_id)

        # ignore
        '''
        newcmp = mpl.colors.LinearSegmentedColormap.from_list('custom blue', 
                                             [(0,         '#0000ff'),
                                              (norm(0.5), '#dddddd'),
                                              (norm(2.0), '#00ffff'),
                                              (1,         '#ff0000')], N=256)
        '''

        norm_grav = mpl.colors.LogNorm( vmin=0.01,vmax=100)
        ok = ~np.isnan(obj.RVg)
        AXX=ax2
        plot=AXX.pcolormesh( XXX,YYY, np.abs(1/obj.RVg),  norm=norm_grav, shading='nearest', cmap=cmap)
        fig.colorbar(plot,label='EK/EG', ax=AXX )
        ax0.set(yscale='log', ylabel=r'$R~~[\rm{AU}]$', xticks=[])
        ax2.set(yscale='log',ylabel=r'$R~~[\rm{AU}]$',xlabel=r'$t/t_{\rm{ff}}$', xlim=ax0.get_xlim())#,ylabel='R[AU]')
        ax1.set(yscale='log',ylabel=r'$R~~[\rm{AU}]$', xticks=[], xlim=ax0.get_xlim())
        if tall==False:
            ax0.set(xlabel='t/tff')
            ax1.set(xlabel='t/tff')
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
        print("Saving!")
        fig.savefig('plots_to_sort/mass_flux_color_%s_c%04d.pdf'%(obj.sim_name, core_id))



# OMG THERE COULD BE DEATH OF FETUS STARS!!!! or simply a merge or accretion...sure.
import track_loader as TL
sim_list = ['m0232']#, 'm0240', 'm0250', 'm0260', 'm0270', 'm0280']
TL.load_tracks(sim_list)

if 'things' not in dir():
    things={}

# use
if 1:
    #new plots.
    for sim in sim_list:
        #core_list = TL.loops[sim].core_by_mode['Alone']
        all_cores=np.unique(TL.loops[sim].tr.core_ids)
        core_list=list(all_cores) 
        print("The core count of target frame %s is %d"%(sim, len(core_list)))
        print("The newborn tags are: ") 
        print(core_list)

        for core_id in core_list:
            if core_id in things:
                continue
            # maybe re-define
            thing = Relaxor(TL.loops[sim])
            thing.run(core_list=[core_id], do_plots=True, frame_list=None)#, tsing=tsing_tool[sim])
            things[core_id]=thing
            plotmancer(thing, tall=True)  #now changed


    if 0:
        plotsb(things, tsing_tool[sim])
    if 0:
        #really, don't overwrite things. please.
        thing_saver(things,sim)
    if 0:
        if 'accreter' not in dir() or True:
            accreter=accretion(things, tsing_tool[sim])
    if 0:
        plotsa(accreter,tsing_tool[sim])
    if 0:
        for core_id in things:
            plotmancer(things[core_id], tsing_tool[sim])
    if 0:
        for core_id in things:
            #plotmancer(things[core_id], tsing_tool[sim])
            calculus(things[core_id],tsing_tool[sim])
if 0:
    plotmancer(thing,tsing_tool[sim])
if 0:
    accretion(thing,tsing_tool[sim])


# NOTES
# use Relaxor, things, and mancer
