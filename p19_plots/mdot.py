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
def calculus(obj, tsing):
        core_id=obj.core_id

        fig,axes=plt.subplots(2,2,figsize=(12,12))
        ax0=axes[0][0]; ax1=axes[0][1]; ax2=axes[1][0]; ax3=axes[1][1]
        #ax0=axes[0][1]; ax1=axes[0][0]
        #ax2=axes[1][0]; ax3=axes[1][1]
        #ax0=axes[1];ax1=axes[0];ax2=axes[2]
        #fig.subplots_adjust(wspace=0)
        rcen = 0.5*(obj.r_bins[1:]+obj.r_bins[:-1])

        cmap=copy.copy(mpl.cm.get_cmap("seismic"))
        cmap.set_bad([0.9]*3)
        XXX,YYY = np.meshgrid( obj.times.flatten(),rcen)

        if 1:
            thalf = 0.5*(obj.times[1:]+obj.times[:-1])
            dt    =     (obj.times[1:]-obj.times[:-1])
            dt.shape=1,dt.size
            MMM = obj.MMM
            dMass = MMM[:,1:]-MMM[:,:-1]
            Mcen  = 0.5*(MMM[:,1:]+MMM[:,:-1])
            dMdt = dMass/dt/Mcen*colors.tff
            ok = ~np.isnan(dMdt)
            maxmax=np.abs(dMdt[ok]).max()
            maxmax=1
            norm2 = mpl.colors.Normalize(-maxmax,maxmax)
            X2,Y2=np.meshgrid(thalf.flatten(),rcen)
            plot0=ax0.pcolormesh(X2 ,Y2, dMdt,  norm=norm2, shading='nearest', cmap=cmap)
            ax0.set(yscale='log',xlabel='t/tff')
            fig.colorbar(plot0,label='From Mass', ax=ax0 )

            ok = ~np.isnan(obj.dMdTf)
            ext = extents(np.abs(obj.dMdTf[ok]))
            norm_flux = mpl.colors.Normalize(vmin=-ext.minmax[1],vmax=ext.minmax[1])
            plot=ax1.pcolormesh(XXX,YYY,obj.dMdTf, norm=norm_flux,shading='nearest',cmap=cmap)
            fig.colorbar(plot, label=r'From Vr',ax=ax1)
            ax1.set(yscale='log', xlabel='t/tff')


        if 0:
            norm_velocity = mpl.colors.Normalize(vmin=-1000,vmax=1000)
            plot=ax0.pcolormesh( XXX,YYY, obj.divv,  norm=norm_velocity, shading='nearest', cmap=cmap)
            plot2=ax1.pcolormesh( XXX,YYY, obj.meanvr,  norm=norm_velocity, shading='nearest', cmap=cmap)
        if 1:
            ok = ~np.isnan(obj.dMdT3d)
            ext = extents(np.abs(obj.dMdT3d[ok]))
            norm_flux = mpl.colors.Normalize(vmin=-ext.minmax[1],vmax=ext.minmax[1])
            plot2=ax2.pcolormesh(XXX,YYY,obj.dMdT3d, norm=norm_flux,shading='nearest',cmap=cmap)
            fig.colorbar(plot2, label=r'vr',ax=ax2)
            ax2.set(yscale='log', xlabel='t/tff', ylabel='R [AU]')


        fig.tight_layout()
        fig.savefig('plots_to_sort/velocity_calculus.png')

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
            mask=ok ;print('kludge mask')
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
        core_list=core_list[4:5]

        core_list = [114]
        for core_id in core_list:
            if core_id in things:
                continue
            thing = Relaxor( TL.loops[sim])
            thing.run(core_list=core_list, do_plots=True, frame_list=None, tsing=tsing_tool[sim])#, r_inflection=anne.inflection[sim])
            things[core_id]=thing

        #plotmancer(thing, tsing_tool[sim])

    if 1:
        for core_id in things:
            #plotmancer(things[core_id], tsing_tool[sim])
            calculus(things[core_id],tsing_tool[sim])

if 0:
    plotmancer(thing,tsing_tool[sim])
