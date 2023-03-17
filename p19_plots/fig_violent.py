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
class Relaxor():
    def __init__(self,this_looper):
        self.this_looper=this_looper
    def run(self,core_list=None, do_plots=True, r_inflection=None, frame_list=None, tsing=None):
        this_looper=self.this_looper

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
            self.RV = np.zeros( [len(r_bins)-1, len(frame_list)])
            self.RVt = np.zeros( [len(r_bins)-1, len(frame_list)])
            self.RVg = np.zeros( [len(r_bins)-1, len(frame_list)])
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

                ok = ~np.isnan(mean_flux)
                self.RV[ok,nframe]=mean_flux[ok]/mass_quant[ok]*colors.tff
                self.RV[~ok,nframe]=np.nan

                ok = ~np.isnan(mass_quant)
                self.RVt[ok,nframe]=mass_quant[ok]
                self.RVt[~ok,nframe]=np.nan

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
def plotmancer(obj, tsing):
        core_id=obj.core_id
        #ax6.set(xlabel='R',xscale='log',ylabel='|V_r|',yscale='log')
        #ax6.set_title('Slope at tend: %0.2f'%pfit[0])
        #ax6.plot( r_bins[1:], 4*np.pi/3*r_bins[1:]**3,'k')
        RV = obj.RV+0
        RVt = obj.RVt+0
        ok = ~np.isnan(RV)
        ext = extents(RV[ok])
        #okt = ~np.isnan(RVt)
        #extt = extents(RVt[ok])
        #for i in range(3):
        #    ax6[i][0].set(ylim=ext.minmax)
        #    ax6[i][1].set(ylim=extt.minmax)
        #fig6.savefig('plots_to_sort/velocity_vs_radius_%s_c%04d.png'%(this_looper.sim_name, core_id))



        fig,axes=plt.subplots(1,3,figsize=(12,3))
        #ax0=axes[0][0]; ax1=axes[0][1]
        #ax2=axes[1][0]; ax3=axes[1][1]
        ax0=axes[1];ax1=axes[0];ax2=axes[2]
        #fig.subplots_adjust(wspace=0)
        rcen = 0.5*(obj.r_bins[1:]+obj.r_bins[:-1])

        #
        # UNITS.  THERES BETTER PLACES FOR THIS.
        #
        #4.9 pc in AU
        rcen *= colors.length_units_au
        #5900 solar masses
        RVt *= colors.mass_units_msun

        XXX,YYY = np.meshgrid( obj.times.flatten(),rcen)
        ax=axes
        #norm = mpl.colors.LogNorm( RV[RV>0].min(), RV.mean())
        ok = ~np.isnan(RV)
        maxmax = np.abs(RV[ok]).max()
        norm = mpl.colors.Normalize( -maxmax,maxmax)
        #pdb.set_trace()
        cmap=copy.copy(mpl.cm.get_cmap("seismic"))
        cmap.set_bad([0.9]*3)

        ok = ~np.isnan(RVt)
        radius_to_cut_off = 1000
        index = np.where( rcen < radius_to_cut_off)[0].max()
        mass_at_the_end = RVt[index,-1]
        fiducial=mass_at_the_end
        Max=RVt[ok].max()
        minner_chicken_dinner = fiducial**2/Max
        norm_mass = mpl.colors.LogNorm(vmin=minner_chicken_dinner,vmax=Max)

        plot=ax0.pcolormesh( XXX,YYY, RV,  norm=norm, shading='nearest', cmap=cmap)
        plot2=ax1.pcolormesh( XXX,YYY, RVt,  norm=norm_mass, shading='nearest', cmap=cmap)
        fig.colorbar(plot, label=r'$\frac{\dot{M}}{M/ t_{ff}}$',ax=ax0)
        fig.colorbar(plot2, label=r'$M(<r) [M_\odot]$',ax=ax1)
        ax0.set(yscale='log', xlabel='t/tff')
        ax1.set(yscale='log', xlabel='t/tff', ylabel='R [AU]')



        newcmp = mpl.colors.LinearSegmentedColormap.from_list('custom blue', 
                                             [(0,         '#0000ff'),
                                              (norm(0.5), '#dddddd'),
                                              (norm(2.0), '#00ffff'),
                                              (1,         '#ff0000')], N=256)

        norm_grav = mpl.colors.LogNorm( vmin=0.01,vmax=100)
        ok = ~np.isnan(obj.RVg)
        print(obj.RVg[ok])
        plot=ax2.pcolormesh( XXX,YYY, np.abs(obj.RVg),  norm=norm_grav, shading='nearest', cmap=cmap)
        fig.colorbar(plot,label='EG/EK', ax=ax2 )
        ax2.set(yscale='log',xlabel='t/tff')#,ylabel='R[AU]')
        #pdb.set_trace()
        
        if 0:
            rrrr = msR.transpose()
            rrrr = rrrr[mask,:]

            ax.plot(times , rrrr, c=[0.5]*3, linewidth=0.1, alpha=0.5)
        if tsing:
            for aaa in axes.flatten():
                aaa.axvline(  tsing.tsing_core[core_id],c='k')
                #ax1.axvline( tsing.tsing_core[core_id],c='k')
                aaa.axvline(  tsing.tend_core[core_id],c='k')
                #ax1.axvline( tsing.tend_core[core_id],c='k')


        fig.tight_layout()
        fig.savefig('plots_to_sort/mass_flux_color_%s_c%04d.pdf'%(obj.this_looper.sim_name, core_id))



if 0:
    import three_loopers_six as TL
    sim_list=['u601','u602','u603']
    sim_list=['u602']
import three_loopers_u500 as TL
sim_list=['u501','u502','u503']
sim_list=['u502']

import tsing
reload(tsing)
if 'tsing_tool' not in dir() or True:
    tsing_tool={}
    for ns,sim in enumerate(sim_list):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()

#import anne
#reload(anne)
##anne.make_inflection()
if 'thing' not in dir():
    for sim in sim_list:
        all_cores=np.unique( TL.loops[sim].tr.core_ids)
        core_list=list(all_cores)
        core_list=None#[323]
        core_list=[323]
        core_list=[25]
        core_list=[74]
        core_list=[114]
        #core_list = TL.loops[sim].core_by_mode['Alone']
        #core_list=core_list[2:3]

        #core_list = [114]
        #core_list=[361]
        #core_list=[8]
        #core_list=[381]
        #core_list=[323]

        thing = Relaxor( TL.loops[sim])
        thing.run(core_list=core_list, do_plots=True, frame_list=None, tsing=tsing_tool[sim])#, r_inflection=anne.inflection[sim])

plotmancer(thing, tsing_tool[sim])
