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
def vel_color(this_looper,core_list=None, do_plots=True, r_inflection=None, frame_list=None, tsing=None):

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
        frame_list=thtr.frames[mask]

    rm = rainbow_map(len(frame_list))
    for core_id in core_list:
        ms = trackage.mini_scrubber(this_looper.tr,core_id,do_velocity=True)
        
        r_bins = np.geomspace( 2e-4, 32/128, 32)
        r_cen = 0.5*(r_bins[1:]+r_bins[:-1])
        RV = np.zeros( [len(r_bins)-1, len(frame_list)])
        RVt = np.zeros( [len(r_bins)-1, len(frame_list)])
        frame_rmap=rainbow_map(len(frame_list))
        fig6,ax6=plt.subplots(3,2)
        ax_index=0
        got_tsing=False
        got_tend=False
        for nframe,frame in enumerate(frame_list):
            print("MassIn on %s c%04d n%04d"%(this_looper.sim_name, core_id, frame))
            ds = this_looper.load(frame)
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


            if 0:
                R1 = sp['radius']
                order = np.argsort(R1)
                vel = []
                for axis in 'xyz':
                    vel.append( sp['velocity_%s'%axis][order][:10].mean())
                scrub = other_scrubber.scrubber(sp, reference_velocity = vel)

                #Get data arrays

                dv = scrub.cell_volume
                RR = scrub.r
                DD = scrub.density
                #dv = np.abs(sp[YT_cell_volume])
                #RR =sp[YT_radius]
                #DD = sp[YT_density]
                vr = scrub.vr_rel
                vt = scrub.vt_rel
            else:
                dv = ms.cell_volume[:,nf]
                RR = ms.r[:,nf]
                DD = ms.density[:,nf]
                vr = ms.vr_rel[:,nf]
                vt = ms.vt2_rel[:,nf]**0.5

            ORDER = np.argsort( RR)
            RR_cuml = RR[ORDER]
            M_local = DD[ORDER]*dv[ORDER]
            M_cuml = np.cumsum( DD[ORDER]*dv[ORDER])
            V_cuml = np.cumsum( dv[ORDER])
            V_local = dv[ORDER]

            #vr_cuml = np.cumsum( scrub.vr_rel[ORDER]*dv[ORDER])
            vr_cumsum = np.cumsum( vr[ORDER]*dv[ORDER])/np.cumsum(dv[ORDER])
            vt_cumsum = np.cumsum( vt[ORDER]*dv[ORDER])/np.cumsum(dv[ORDER])


            digitized = np.digitize( RR_cuml, r_bins)
            #vr_mean  =nar([ (np.abs(M_cuml[ digitized == i])*dv[digitized==i]).sum()/(dv[digitized==i].sum()) if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])
            #vr_mean  =nar([ M_cuml[ digitized == i].max() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])
            vr_mean  =nar([ vr_cumsum[ digitized == i].mean() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])
            vt_mean  =nar([ vt_cumsum[ digitized == i].mean() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])
            ok = ~np.isnan(vr_mean)
            RV[ok,nframe]=vr_mean[ok]
            RV[~ok,nframe]=np.nan
            RVt[ok,nframe]=vt_mean[ok]
            RVt[~ok,nframe]=np.nan
            r_bins_c = 0.5*(r_bins[1:]+r_bins[:-1])
            #ax6.plot(r_bins_c, vr_mean, c = frame_rmap(nframe))
            linewidth=0.3
            c=frame_rmap(nframe)
            if 0:
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


                ##ax6.plot(RR_cuml, gaussian_filter(M_cuml,1), c=c,linewidth=linewidth)
                ##ax6.plot(RR_cuml[1:], M_cuml[1:]/RR_cuml[1:]**3, c=c,linewidth=linewidth)
                ax6[ax_index][0].plot( r_cen[ok], vr_mean[ok], c=c, linewidth=linewidth)
                ax6[ax_index][1].plot( r_cen[ok], vt_mean[ok], c=c, linewidth=linewidth)
                #ax6[0].hist( RR_cuml.v, histtype='step', color=c)
        #ax6.set(xlabel='R',xscale='log',ylabel='|V_r|',yscale='log')
        #ax6.set_title('Slope at tend: %0.2f'%pfit[0])
        #ax6.plot( r_bins[1:], 4*np.pi/3*r_bins[1:]**3,'k')
        ok = ~np.isnan(RV)
        ext = extents(RV[ok])
        okt = ~np.isnan(RVt)
        extt = extents(RVt[ok])
        for i in range(3):
            ax6[i][0].set(ylim=ext.minmax)
            ax6[i][1].set(ylim=extt.minmax)
        fig6.savefig('plots_to_sort/velocity_vs_radius_%s_c%04d.png'%(this_looper.sim_name, core_id))



        fig,axes=plt.subplots(1,2,figsize=(12,8))
        #fig.subplots_adjust(wspace=0)
        rcen = 0.5*(r_bins[1:]+r_bins[:-1])
        XXX,YYY = np.meshgrid( times.flatten(),rcen)
        ax=axes
        #norm = mpl.colors.LogNorm( RV[RV>0].min(), RV.mean())
        ok = ~np.isnan(RV)
        maxmax = np.abs(RV[ok]).max()
        maxmax=max([maxmax, RVt[ok].max()])
        norm = mpl.colors.Normalize( -maxmax,maxmax)
        cmap=copy.copy(mpl.cm.get_cmap("seismic"))
        cmap.set_bad([0.9]*3)

        plot=ax[0].pcolormesh( XXX,YYY, RV,  norm=norm, shading='nearest', cmap=cmap)
        plot2=ax[1].pcolormesh( XXX,YYY, RVt,  norm=norm, shading='nearest', cmap=cmap)
        fig.colorbar(plot, label='<V_r>(<r)',ax=ax[0])
        fig.colorbar(plot2, label='<V_t>(<r)',ax=ax[1])
        ax[0].set(yscale='log', xlabel='t/tff', ylabel='R [code units]')
        ax[1].set(yscale='log', xlabel='t/tff', ylabel='R [code units]')
        
        if 0:
            rrrr = msR.transpose()
            rrrr = rrrr[mask,:]

            ax.plot(times , rrrr, c=[0.5]*3, linewidth=0.1, alpha=0.5)
        if tsing:
            for aaa in ax:
                aaa.axvline(  tsing.tsing_core[core_id],c='k')
                #ax1.axvline( tsing.tsing_core[core_id],c='k')
                aaa.axvline(  tsing.tend_core[core_id],c='k')
                #ax1.axvline( tsing.tend_core[core_id],c='k')


        fig.savefig('plots_to_sort/velocity_color_%s_c%04d.png'%(this_looper.sim_name, core_id))



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
if 'RV' not in dir() or True:
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
        TF=TL.loops[sim].target_frame
        all_frames = TL.loops[sim].tr.frames
        nframes = all_frames.size
        frame_list = all_frames #[0, TF, all_frames[int(nframes/2)]]
        #frame_list = [100]

        RV=vel_color( TL.loops[sim],core_list=core_list, do_plots=True, frame_list=None, tsing=tsing_tool[sim])#, r_inflection=anne.inflection[sim])
