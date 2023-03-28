from starter2 import *
import xtra_energy

from scipy.optimize import curve_fit
from scipy import stats
from scipy.ndimage import gaussian_filter
import core_proj_three
reload(core_proj_three)
import other_scrubber
reload(other_scrubber)
#import three_loopers_six as TL
import camera_path
import three_loopers_u500 as TL
sim_list=['u501','u502','u503']
sim_list=['u502']
plt.close('all')
if 0:
    for sim in sim_list:
        loop = TL.loops[sim]
        camera = camera_path.camera_1( loop, 'smooth_zoom_2')
        core_proj_three.core_proj_multiple(loop,axis_list=[0],core_list=[74],frame_list=[0,10],camera=camera, main_core=74)


class multi_profile():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
        self.profiles={}

    def run(self, core_list=None, frame_list=None, tsing=None):
        this_looper=self.this_looper
        thtr=this_looper.tr
        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)
        def get_time_index(time):
            index=np.argmin( np.abs( thtr.times/colors.tff-time))
            return index
        Nplots = 4
        Ntimes = 4
        fig,axes = plt.subplots(Nplots,Ntimes, figsize=(8,12))
        fig.subplots_adjust(hspace=0,wspace=0)
        ext = [extents() for n in range(Nplots+1)]
        if Nplots > 4:
            rho_ext = extents()
            for core_id in core_list:
                ms = trackage.mini_scrubber(this_looper.tr,core_id)
                rho_ext(ms.density)
        for core_id in core_list:
            self.profiles[core_id]={}
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            self.cores_used.append(core_id)

            frame_mask = np.zeros_like(thtr.times, dtype='bool')
            frame_mask[0]=True
            frame_mask[get_time_index(0.9*tsing.tsing_core[core_id])]=True
            frame_mask[get_time_index(tsing.tsing_core[core_id])]=True
            frame_mask[get_time_index(tsing.tend_core[core_id])]=True
            frame_list=thtr.frames[frame_mask]
            rm = rainbow_map(len(frame_list))
            axes[0][0].set_title(r'$t=0$')
            axes[0][1].set_title(r'$t=0.9 t_{\rm{sing}}$')
            axes[0][2].set_title(r'$t=t_{\rm{sing}}$')
            axes[0][3].set_title(r'$t=t_{\rm{end}}$')

            img_collector=[]
            for nframe,frame in enumerate(frame_list):
                self.profiles[core_id][frame]={}
                print("profile on %s c%04d n%04d"%(this_looper.sim_name, core_id, frame))
                ds = this_looper.load(frame)
                xtra_energy.add_energies(ds)
                xtra_energy.add_gdotgradrho(ds)
                nf = np.where( this_looper.tr.frames == frame)[0][0]

                center = nar([ms.mean_xc[nf], ms.mean_yc[nf],ms.mean_zc[nf]])
                msR = ms.rc
                msR[ msR<1/2048]=1/2048
                
                MaxRadius=msR[:,nf].max()
                Radius = max([8.0/128, MaxRadius])
                rsph = ds.arr(Radius,'code_length')
                sph = ds.sphere(center,rsph)

                dv = sph[YT_cell_volume]
                RR = sph['radius']
                DD = sph[YT_density]
                ORDER = np.argsort( RR)
                rho_sort = DD[ORDER]
                RR_sort = RR[ORDER]
                dv_sort = dv[ORDER]
                M_cuml = np.cumsum( DD[ORDER]*dv[ORDER])
                d2_cuml = np.cumsum( DD[ORDER]**2*dv[ORDER])
                V_cuml = np.cumsum( dv[ORDER])
                if 1:
                    vel = []
                    for axis in 'xyz':
                        vel.append( sph['velocity_%s'%axis][ORDER][:10].mean())
                        #vel.append( sp['velocity_%s'%axis][ORDER].mean())
                        #vel.append( (rho_sort*sp['velocity_%s'%axis][ORDER]*dv_sort)[:30].sum()/M_cuml)
                        #vel.append( (rho_sort*sph['velocity_%s'%axis][ORDER]*dv_sort)[:30].sum()/(rho_sort*dv_sort)[:30].sum())
                        #vel.append(0)
                    scrub = other_scrubber.scrubber(sph, reference_velocity = vel)
                    #scrub0 = other_scrubber.scrubber(sp, reference_velocity = [0.0]*3)
                    scrub.compute_ke_rel()
                EG = sph[YT_grav_energy_2]
                EK = scrub.ke_rel
                EG_cuml = np.abs(np.cumsum( EG[ORDER]*dv[ORDER]))/V_cuml
                EK_cuml = np.cumsum( EK[ORDER]*dv[ORDER])/V_cuml


                #pw = yt.ProjectionPlot(ds, ax, YT_density, data_source=sph, center=center, origin='window', weight_field=weight_field)
                #img_collector.append( pw.data_source.to_frb( 2*sph.radius,512))
                #axes[0][nframe].imshow( img_collector[-1][YT_density], norm=mpl.colors.LogNorm())
                #axes[0][nframe].set(xticks=[], yticks=[])

                args = {'linewidth':0.2, 'c':[0.5]*4}
                ext[-1](RR_sort)
                #axes[0][nframe].plot( RR_sort, rho_sort, c=[0.5]*4)
                axes[0][nframe].plot( RR_sort, colors.density_units* M_cuml/V_cuml, **args)
                ext[0](colors.density_units*M_cuml/V_cuml)
                self.profiles[core_id][frame]['rho']=colors.density_units* M_cuml/V_cuml


                vr = scrub.vr_rel
                vr_cumsum = np.cumsum( vr[ORDER]*dv_sort)/V_cuml
                #vr_cumsum = np.abs(np.cumsum( vr[ORDER]*dv_sort)/V_cuml)
                #vr_cumsum = np.abs(np.cumsum( rho_sort*vr[ORDER]*dv_sort)/V_cuml)
                self.profiles[core_id][frame]['R']=RR_sort
                self.profiles[core_id][frame]['vr_cumsum']=vr_cumsum
                axes[1][nframe].plot(RR_sort, vr_cumsum, **args)
                ext[1](vr_cumsum)

                vt = scrub.vt_rel
                vt_cumsum = np.cumsum( vt[ORDER]*dv_sort)/V_cuml
                self.profiles[core_id][frame]['vt_cumsum']=vt_cumsum
                axes[2][nframe].plot(RR_sort, vt_cumsum, **args)
                ext[2](vt_cumsum)
                #if 'did_one_plot' not in dir():
                #    did_one_plot=False
                #if nframe==2 and not did_one_plot:
                #    did_one_plot=True
                #    print("WOOO")
                #    axes[2][nframe].plot(RR_sort, RR_sort*vt_cumsum.max()/RR_sort.max())

                #axes[3][nframe].plot(RR_sort,EG_cuml/EK_cuml,**args)
                axes[3][nframe].plot(RR_sort,EK_cuml/EG_cuml,**args)
                self.profiles[core_id][frame]['energy']=EK_cuml/EG_cuml
                self.profiles[core_id][frame]['EGmean']=EG_cuml
                self.profiles[core_id][frame]['EKmean']=EK_cuml
                ext[3](EK_cuml/EG_cuml)

                if Nplots > 4:
                    bins = np.geomspace(rho_ext.minmax[0], rho_ext.minmax[1],128)
                    rho_to_hist = DD+0
                    rho_to_hist.sort()
                    cuml = np.arange(rho_to_hist.size)/rho_to_hist.size
                    #axes[4][0].plot( rho_to_hist, cuml)
                    axes[4][nframe].hist( DD.v, weights=dv.v, histtype='step', bins=bins, density=True)
                if 'hru' not in dir():
                    hru=0
            #ooo='plots_to_sort/incremental_c%04d'%core_id
            #fig.savefig(ooo)
            #print(ooo)

        for ax in axes[0]:
            ax.set(xscale='log',yscale='log',ylim=ext[0].minmax, ylabel=r'$<\rho>(<r)$ [\rm{g/cc}]$', xlim=ext[-1].minmax)
        for ax in axes[1]:
            ax.set(xscale='log',yscale='linear',ylim=ext[1].minmax, ylabel=r'$<v_r>(<r)$', xlim=ext[-1].minmax)
            ax.axhline(0,c=[0.5]*4)
        for nframe,ax in enumerate(axes[2]):
            ax.set(xscale='log',yscale='linear',ylim=ext[2].minmax, ylabel=r'$<v_t>(<r)$', xlim=ext[-1].minmax)
            if nframe==2:
                RR_sort = np.geomspace(1e-4,1e1)
                #ax.plot(RR_sort, RR_sort*vt_cumsum[-1]/RR_sort[-1], c='r')
                ax.plot(RR_sort, 4*(RR_sort/0.1)**0.5, c='r')

        for nframe,ax in enumerate(axes[3]):
            #THIS ONE
            #ax.set(xscale='log',yscale='log',ylim=ext[3].minmax, ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)

            ax.set(xscale='log',yscale='linear',ylim=[0.5,5], ylabel=r'$Ek(<r)/Eg(<r)$', xlim=ext[-1].minmax)

            #ax.set(xscale='log',yscale='log',ylim=ext[3].minmax, ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)
            #ax.set(xscale='log',yscale='linear',ylim=ext[3].minmax, ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)
            #ax.set(xscale='log',yscale='linear',ylim=[0,3], ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)
            #ax.set(xscale='log',yscale='log',ylim=[0,3], ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)
            #ax.set(xscale='log',yscale='linear',ylim=[0,3], ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)
            ax.axhline(0.5,c=[0.5]*3)
            ax.axhline(2.0,c=[0.5]*3)
            
            if nframe==3 and False:
                RR_sort = np.geomspace(1e-4,1e1)
                rmin=RR_sort.min()
                X = np.log((RR_sort/rmin))
                x_scale = np.log(X)
                x0, x1 = np.log(1e-3/rmin), np.log(1e-2/rmin)
                y0, y1 = 0.5, 2
                m = (y1-y0)/(x1-x0)
                b = y0 - m*x0
                print("M",m)
                print("uM", 1.5/-np.log(10))
                print("b %0.5f"%b, y1-y0, x1-x0)
                they=m*X+b
                axes[3][nframe].plot(RR_sort, they, c='r')


        if Nplots>4:
            for ax in axes[4]:
                ax.set(xscale='log',yscale='log')
        for row in axes:
            for ax in row[1:]:
                ax.set(ylabel='', yticks=[])
        for ax in axes[-1]:
            ax.set(xlabel='r')
        fig.savefig('plots_to_sort/radial_profile_horse_%s'%(this_looper.sim_name))

def geplotter(obj):
    Nplots=3
    Ntimes=4
    fig,axes = plt.subplots(Nplots,Ntimes, figsize=(8,12))
    fig.subplots_adjust(hspace=0,wspace=0)
    axes[0][0].set_title(r'$t=0$')
    axes[0][1].set_title(r'$t=0.9 t_{\rm{sing}}$')
    axes[0][2].set_title(r'$t=t_{\rm{sing}}$')
    axes[0][3].set_title(r'$t=t_{\rm{end}}$')
    row_dict={'EGmean':0,'EKmean':1,'energy':2}
    
    ext = [extents() for n in range(Nplots+1)]
    args = {'linewidth':0.2, 'c':[0.5]*4}
    profs=['EGmean','EKmean','energy']
    #profs=['energy']
    for nprofile, profile in enumerate(profs):
        print(profile)

        core_list=list(obj.profiles.keys())
        #core_list=core_list[:5]
        for core_id in core_list:
            row = row_dict[profile]
            frames = sorted(list(obj.profiles[core_id].keys()))
            for nframe,frame in enumerate(frames):
                if nframe >0:
                    continue
                ax = axes[row][nframe]
                R = nar(obj.profiles[core_id][frame]['R'].v)
                y = nar(obj.profiles[core_id][frame][profile].v)
                if profile == 'EKmean':
                    #y=y/R
                    rbins = np.geomspace(1e-4,R.max(),64)
                    ybinned, edge, num = stats.binned_statistic(R,y,statistic='mean',bins=rbins)
                    ok=~np.isnan(ybinned)
                    ybinned = gaussian_filter(ybinned[ok],3)
                    rcen=0.5*(rbins[1:]+rbins[:-1])
                    pfit=np.polyfit(np.log10(rcen[ok]),np.log10(ybinned),1)
                    ax.plot(rcen,10**(pfit[0]*np.log10(rcen)+pfit[1]),c=[0.5]*3,linewidth=0.1,linestyle='--')
                    ax.plot(rcen[ok],ybinned,**args)
                    print(pfit)
                else:
                    ax.plot(R,y,**args)
                ext[-1](R)
                ext[row](y)
    for ax in axes[0]:
        ax.set(xscale='log',yscale='log',ylim=ext[0].minmax, ylabel=r'$EG$', xlim=ext[-1].minmax)
    for ax in axes[1]:
        ax.set(xscale='log',yscale='log',ylim=ext[1].minmax, ylabel=r'$EK$', xlim=ext[-1].minmax)
    for ax in axes[2]:
        ax.set(xscale='log',yscale='log',ylim=ext[2].minmax, ylabel=r'$EK/EG$', xlim=ext[-1].minmax)
    for row in axes:
        for ax in row[1:]:
            ax.set(ylabel='', yticks=[])
    for ax in axes[-1]:
        ax.set(xlabel='r')
    fig.savefig('plots_to_sort/energies.png')

def replotter(obj):
    Nplots=4
    Ntimes=4
    fig,axes = plt.subplots(Nplots,Ntimes, figsize=(8,12))
    fig.subplots_adjust(hspace=0,wspace=0)
    axes[0][0].set_title(r'$t=0$')
    axes[0][1].set_title(r'$t=0.9 t_{\rm{sing}}$')
    axes[0][2].set_title(r'$t=t_{\rm{sing}}$')
    axes[0][3].set_title(r'$t=t_{\rm{end}}$')
    row_dict={'rho':0,'vr_cumsum':1,'vt_cumsum':2,'energy':3}
    
    ext = [extents() for n in range(Nplots+1)]
    args = {'linewidth':0.2, 'c':[0.5]*4}
    profs=['rho','vr_cumsum','vt_cumsum','energy']
    #profs=['energy']
    fit_a=[]
    fit_p=[]
    for nprofile, profile in enumerate(profs):
        print(profile)

        core_list=list(obj.profiles.keys())
        #core_list=core_list[:5]
        for core_id in core_list:
            row = row_dict[profile]
            frames = sorted(list(obj.profiles[core_id].keys()))
            for nframe,frame in enumerate(frames):
                ax = axes[row][nframe]
                R = nar(obj.profiles[core_id][frame]['R'].v)
                ext[-1](R)
                y = nar(obj.profiles[core_id][frame][profile].v)

                if profile != 'energy':
                    y=y
                    ext[row](y)
                    ax.plot( R, y, **args)
                if profile=='energy':
                    #def fun(R,A,P):
                    #    return nar(A/(R/1e-3)**2)+P
                    y=1/y
                    if nframe>0:
                        ax.plot( R, y, **args)
                    if nframe==0:
                        ax.plot( R, y, **args)
                        if 0:
                            rbins = np.geomspace(1e-4,R.max(),64)
                            ybinned, edge, num = stats.binned_statistic(R,y,statistic='mean',bins=rbins)
                            ybinned = gaussian_filter(ybinned,3)

                            ok = ~np.isnan(ybinned)
                            rcen = 0.5*(rbins[1:]+rbins[:-1])
                            ax.plot(rcen[ok],ybinned[ok],**args)
                            pfit=np.polyfit(np.log10(rcen[ok]), np.log10(ybinned[ok]),1)
                            print(pfit)
                        if 0:
                            ax.plot(R, (R/0.1)**(-2))

                    ext[row](y)
                    if nframe==3:
                        rmin=1e-4
                        def fun(R,A,P):
                            return A*np.log(R/rmin)+P
                        ok = (R>1e-3)*(R<1e-2)
                        rbins = np.geomspace(5e-4,1e-2)
                        ybinned, edge, num = stats.binned_statistic(R[ok],y[ok],statistic='mean',bins=rbins)
                        ok = ~np.isnan(ybinned)
                        rcen = 0.5*(rbins[1:]+rbins[:-1])

                        popt, pcov = curve_fit(fun,rcen[ok],ybinned[ok],p0=[2,-1])
                        fit_a.append(popt[0])
                        fit_p.append(popt[1])
                        #print(y)
                        #print(popt)
                        #ax.plot(R, fun(R,*popt),c='r')
                        #ax.plot( rcen, ybinned, c='k')
                        #pfit = np.polyfit(np.log(R[ok]/R[ok].min()), y[ok],1)
                        #ax.plot(R, np.exp( pfit[0]*np.log(R[ok]/R[ok].min())+pfit[1]))
                        #print(pfit)

                        #fig2,ax2=plt.subplots(1,1)
                        #ax2.plot(sorted(R))
                        #ax2.plot(sorted(y)[4000:])
                        #pdb.set_trace()
                        #print(sorted(y)[:10])
                        #fig2.savefig('plots_to_sort/argus.png')

                        if 0:
                            RR_sort = np.geomspace(1e-4,1e1)
                            rmin=RR_sort.min()
                            X = np.log((RR_sort/rmin))
                            x_scale = np.log(X)
                            x0, x1 = np.log(1e-3/rmin), np.log(1e-2/rmin)
                            y0, y1 = 0.5, 2
                            m = (y1-y0)/(x1-x0)
                            b = y0 - m*x0
                            print("M",m, "M*log10", m*np.log(10), "rmin",rmin)
                            print("uM", 1.5/-np.log(10))
                            print("b %0.5f"%b, y1-y0, x1-x0)
                            they=m*X+b
                        #ax.plot(RR_sort, they, c='r')
    if 1:
        fig2,ax2=plt.subplots(1,2)
        a0=ax2[0];a1=ax2[1]
        a0.hist(fit_a)
        a1.hist(fit_p)
        fig2.savefig('plots_to_sort/a_p.png')



    print('set')
    for nax,ax in enumerate(axes[0]):
        ax.set(xscale='log',yscale='log',ylim=ext[0].minmax, ylabel=r'fff$<\rho>(<r)$ [\rm{g/cc}]$', xlim=ext[-1].minmax)
        if nax==3:
            r = np.geomspace(1e-2,0.2)
            ax.plot( r,1e6*(r/0.1)**-2,c='r')
        
    for ax in axes[1]:
        ax.set(xscale='log',yscale='linear',ylim=ext[1].minmax, ylabel=r'$<v_r>(<r)$', xlim=ext[-1].minmax)
        ax.axhline(0,c=[0.5]*4)

        if nax==2:
            r=np.geomspace(1e-3,1e-2)
            y = -1*np.log(r/r.min())
    for nframe,ax in enumerate(axes[2]):
        ax.set(xscale='log',yscale='linear',ylim=ext[2].minmax, ylabel=r'$<v_t>(<r)$', xlim=ext[-1].minmax)
        if nframe==2:
            RR_sort = np.geomspace(1e-4,1e1)
            #ax.plot(RR_sort, RR_sort*vt_cumsum[-1]/RR_sort[-1], c='r')
            ax.plot(RR_sort, 4*(RR_sort/0.1)**0.5, c='r')

    for nframe,ax in enumerate(axes[3]):
        #THIS ONE
        ax.set(xscale='log',yscale='log',ylim=ext[3].minmax, ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)

        #ax.set(xscale='log',yscale='linear',ylim=[0.1,5], ylabel=r'$Ek(<r)/Eg(<r)$', xlim=ext[-1].minmax)
        #ax.plot(R, nar(2/(R/1e-3)**2)+0.5)
        #ax.plot(R, -np.log(R/R.max()))

        #ax.set(xscale='log',yscale='log',ylim=ext[3].minmax, ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)
        #ax.set(xscale='log',yscale='linear',ylim=ext[3].minmax, ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)
        #ax.set(xscale='log',yscale='linear',ylim=[0,3], ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)
        #ax.set(xscale='log',yscale='log',ylim=[0,3], ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)
        #ax.set(xscale='log',yscale='linear',ylim=[0,3], ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)
        ax.axhline(0.5,c=[0.5]*3)
        ax.axhline(2.0,c=[0.5]*3)

        if nframe==3:
            RR_sort = np.geomspace(1e-4,1e1)
            mean_a = nar(fit_a).mean()
            mean_p = nar(fit_p).mean()
            ax.plot(rbins, mean_a*np.log(rbins/1e-4)+mean_p, c='r')
        if nframe==0:
            r2 = np.geomspace(0.05,0.5,32)
            y = 10*(r2/0.1)**-1
            ax.plot(r2,y,c='r')
        
        if nframe==3 and False:
            RR_sort = np.geomspace(1e-4,1e1)
            rmin=RR_sort.min()
            X = np.log((RR_sort/rmin))
            x_scale = np.log(X)
            x0, x1 = np.log(1e-3/rmin), np.log(1e-2/rmin)
            y0, y1 = 0.5, 2
            m = (y1-y0)/(x1-x0)
            b = y0 - m*x0
            print("M",m)
            print("uM", 1.5/-np.log(10))
            print("b %0.5f"%b, y1-y0, x1-x0)
            they=m*X+b
            axes[3][nframe].plot(RR_sort, they, c='r')


    if Nplots>4:
        for ax in axes[4]:
            ax.set(xscale='log',yscale='log')
    for row in axes:
        for ax in row[1:]:
            ax.set(ylabel='', yticks=[])
    for ax in axes[-1]:
        ax.set(xlabel='r')
    print('saving')
    fig.savefig('plots_to_sort/radial_profile_horse_%s'%(obj.this_looper.sim_name))




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
if 'mp' not in dir():
    for sim in sim_list:
        all_cores=np.unique( TL.loops[sim].tr.core_ids)
        core_list=list(all_cores)
        core_list=None#[323]
        core_list=[323]
        core_list=[25]
        core_list=[74]
        core_list = TL.loops[sim].core_by_mode['Alone']
        #core_list=core_list[:2]
        #core_list=core_list[:10]
        #core_list=[68]
        #core_list
        #core_list=core_list[10:]
        #core_list=[114]


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

        mp=multi_profile(TL.loops[sim])
        mp.run(core_list=core_list,tsing=tsing_tool[sim])#, r_inflection=anne.inflection[sim])
if 1:
    replotter(mp)
if 0:
    geplotter(mp)
