from starter2 import *
import xtra_energy
import three_loopers_u500 as TL
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import pcolormesh_helper as pch
import colors
import movie_frames 
from collections import defaultdict
G = colors.G
class vel_color():
    def __init__(self):
        self.alpha_time=defaultdict(list)
        self.alpha_2_time=defaultdict(list)
        self.cores_used=[]
    def run(self,this_looper,core_list=None, do_plots=True, r_inflection=None, frame_list=None, tsing=None):

        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        thtr=this_looper.tr
        if frame_list is None:
            mask = movie_frames.quantized_mask(this_looper).flatten()
            ok = np.zeros_like(mask)
            ok[::10] = mask[::10]
            #mask=ok ;print('kludge mask')
            times=thtr.times[mask]+0 #the zero makes a copy
            self.times=times+0
            times.shape=times.size,1
            times=times/colors.tff
            self.times=times.flatten()
            frame_list=thtr.frames[mask]
            rm = rainbow_map(len(frame_list))




        #fig5,ax5=plt.subplots(1,1)
        for core_id in core_list:
            if 0:
                just_tsing = np.zeros_like(thtr.times, dtype='bool')
                index=np.argmin( np.abs( thtr.times/colors.tff-tsing.tend_core[core_id]))
                print(index)
                just_tsing[index]=True
                times = thtr.times[just_tsing]+0
                times.shape=times.size,1
                frame_list=thtr.frames[just_tsing]
                rm = rainbow_map(len(frame_list))
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            self.cores_used.append(core_id)
            
            r_bins = np.geomspace( 2e-4, 32/128, 32)
            r_cen = 0.5*(r_bins[1:]+r_bins[:-1])
            RV = np.zeros( [len(r_bins)-1, len(frame_list)])
            frame_rmap=rainbow_map(len(frame_list))
            fig6,ax6=plt.subplots(1,1)
            got_tsing=False
            got_tend=False
            for nframe,frame in enumerate(frame_list):
                print("profile on %s c%04d n%04d"%(this_looper.sim_name, core_id, frame))
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


                R1 = sp['radius']
                order = np.argsort(R1)
                #vel = []
                #for axis in 'xyz':
                #    vel.append( sp['velocity_%s'%axis][order][:10].mean())
                #scrub = other_scrubber.scrubber(sp, reference_velocity = vel)

                #Get data arrays

                dv = sp[YT_cell_volume]
                RR = sp['radius']
                DD = sp[YT_density]
                #dv = np.abs(sp[YT_cell_volume])
                #RR =sp[YT_radius]
                #DD = sp[YT_density]

                ORDER = np.argsort( RR)
                RR_cuml = RR[ORDER]
                M_local = DD[ORDER]*dv[ORDER]
                M_cuml = np.cumsum( DD[ORDER]*dv[ORDER])
                M2_cuml = np.cumsum( DD[ORDER]**2*dv[ORDER])
                V_cuml = np.cumsum( dv[ORDER])
                V_local = dv[ORDER]
                rho_cuml = M_cuml/V_cuml

                #digitized = np.digitize( RR_cuml, r_bins)
                #rho  =nar([ rho_cuml[ digitized == i].max() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])
                #ok = ~np.isnan(rho)
                #RV[ok,nframe]=rho[ok]
                #RV[~ok,nframe]=np.nan
                r_bins_c = 0.5*(r_bins[1:]+r_bins[:-1])
                #ax6.plot(r_bins_c, vr_mean, c = frame_rmap(nframe))
                linewidth=0.3
                c=frame_rmap(nframe)
                #pfit = np.polyfit( np.log10(RR_cuml[1:]), np.log10(M_cuml[1:]/RR_cuml[1:]**3),1)
                pfit = np.polyfit( np.log10(RR_cuml[1:]), np.log10(M_cuml[1:]),1)
                self.alpha_time[core_id].append(pfit[0])
                pfit2 = np.polyfit( np.log10(RR_cuml[1:]), np.log10(M2_cuml[1:]),1)
                self.alpha_2_time[core_id].append(pfit2[0])

                #from scipy.ndimage import gaussian_filter
#
#                ax6.plot( RR_cuml[1:], rho_cuml[1:], c=c, linewidth=linewidth)
#                ax6.plot( RR_cuml[1:], 10**(-2*np.log10(RR_cuml[1:]) + pfit[1]))
                
#            ax6.set(xlabel='R',xscale='log',ylabel='|rho|',yscale='log')
            #ax6.set_title('Slope at tend: %0.2f'%pfit[0])
            #ax6.plot( r_bins[1:], 4*np.pi/3*r_bins[1:]**3,'k')
            #ok = ~np.isnan(RV)
            #ext = extents(RV[ok])
            #for i in range(3):
            #    ax6.set(ylim=ext.minmax)
            #fig6.savefig('plots_to_sort/mean_density_vs_radius_%s_c%04d.png'%(this_looper.sim_name, core_id))

            #ax5.plot(times.flatten(), alpha_time[core_id])
        #fig5.savefig('plots_to_sort/several_alpha_%s'%(this_looper.sim_name))



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
if 'RV' not in dir():
    for sim in sim_list:
        all_cores=np.unique( TL.loops[sim].tr.core_ids)
        core_list=list(all_cores)
        core_list=None#[323]
        core_list=[323]
        core_list=[25]
        core_list=[74]
        core_list = TL.loops[sim].core_by_mode['Alone']
        #core_list=core_list[:10]

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

        RV=vel_color()
        RV.run(TL.loops[sim],core_list=core_list, do_plots=True, frame_list=None, tsing=tsing_tool[sim])#, r_inflection=anne.inflection[sim])
if 0:
    #RV=RVsave
    fig,ax=plt.subplots(1,1)
    for core_id in RV.cores_used:
        ax.plot( RV.times/tsing_tool[sim].tend_core[core_id], RV.alpha_time[core_id])
        ax.axhline(-2)
        ax.axvline(1)
        #print(1-tsing_tool[sim].tend_core[core_id]/tsing_tool[sim].tsing_core[core_id])
    ax.set(xlabel='t/tff',ylabel=r'$\alpha_{\rho,r}$')
    fig.savefig('plots_to_sort/alphas_tend_%s.png'%sim)
    fig,ax=plt.subplots(1,1)
    for core_id in RV.cores_used:
        ax.plot( RV.times/tsing_tool[sim].tend_core[core_id], RV.alpha_2_time[core_id])
        #print(1-tsing_tool[sim].tend_core[core_id]/tsing_tool[sim].tsing_core[core_id])
    ax.set(xlabel='t/tff',ylabel=r'$\alpha_{\rho^2,r}$')
    fig.savefig('plots_to_sort/alphas_2_tend_%s.png'%sim)
if 0:
    alpha_by_time=defaultdict(list)
    alpha2_by_time=defaultdict(list)
    for core_id in RV.cores_used:
        tsing = tsing_tool[sim].tsing_core[core_id]
        tend = tsing_tool[sim].tend_core[core_id]
        for nt,t in enumerate([0.5*tsing, 0.9*tsing, tsing, tend]):
            index = np.argmin( np.abs(RV.times-t))
            alpha2_by_time[nt].append( RV.alpha_2_time[core_id][index])
            alpha_by_time[nt].append( RV.alpha_time[core_id][index])
    fig,ax=plt.subplots(1,2)
    rm = rainbow_map(4)
    for nt in range(4):
        ax[0].hist(alpha_by_time[nt], histtype='step',  color=rm(nt), cumulative=True)
        ax[1].hist(alpha2_by_time[nt], histtype='step', color=rm(nt), cumulative=True)
    ax[0].set(ylabel='N',xlabel='alpha')
    ax[1].set(ylabel='N',xlabel='alpha2')
    fig.savefig('plots_to_sort/alpha_dist_by_tsing.png')
import multiplots
if 1:
    alpha_by_time=defaultdict(list)
    alpha2_by_time=defaultdict(list)
    for core_id in RV.cores_used:
        tsing = tsing_tool[sim].tsing_core[core_id]
        tend = tsing_tool[sim].tend_core[core_id]
        for nt,t in enumerate([0.5*tsing, 0.9*tsing, tsing, tend]):
            index = np.argmin( np.abs(RV.times-t))
            alpha2_by_time[nt].append( RV.alpha_2_time[core_id][index])
            alpha_by_time[nt].append( RV.alpha_time[core_id][index])
    plt.close('all')
    plt.clf()
    fig, ax1,axtop,axright=multiplots.three_way_bean(figsize=(4,4), left=0.15, width=0.62, bottom=0.11, height=0.62, histdepth=0.02)
    rm = rainbow_map(4)
    for nt in range(4):
        ax1.scatter(alpha_by_time[nt], alpha2_by_time[nt], color=rm(nt))
        axtop.hist( alpha_by_time[nt], histtype='step',cumulative=False, color=rm(nt))
        axright.hist( alpha2_by_time[nt], histtype='step',cumulative=False, color=rm(nt), orientation='horizontal')
    ax1.set(ylabel='alpha2',xlabel='alpha')
    fig.savefig('plots_to_sort/alpha1_alpha2.png')

plt.close('all')


