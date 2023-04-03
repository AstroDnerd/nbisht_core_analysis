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

    def run(self, core_list=None, frame_list=None, tsing=None, do_plots=False, timescale=0):
        self.timescale=timescale
        this_looper=self.this_looper
        thtr=this_looper.tr
        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)
        def get_time_index(time):
            index=np.argmin( np.abs( thtr.times/colors.tff-time))
            return index
        Nplots = 4
        Ntimes = 6
        if do_plots:
            fig,axes = plt.subplots(Nplots,Ntimes, figsize=(12,12))
            fig.subplots_adjust(hspace=0,wspace=0)
        ext = [extents() for n in range(Nplots+1)]
        for core_id in core_list:
            self.profiles[core_id]={}
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            self.cores_used.append(core_id)

            frame_mask = np.zeros_like(thtr.times, dtype='bool')
            if self.timescale==1:
                self.titles=[]
                for theta in np.linspace(0,1,Ntimes):
                    t = (1-theta)*tsing.tsing_core[core_id]+theta*tsing.tend_core[core_id]
                    self.titles.append(r'%0.2f'%theta)
                    index=get_time_index(t)
                    frame_mask[index]=True


            if self.timescale==0:
                frame_mask[0]=True
                frame_mask[get_time_index(0.25*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(0.5*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(0.75*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(tsing.tend_core[core_id])]=True
                #half_collapse = 0.5*(tsing.tsing_core[core_id]+tsing.tend_core[core_id])
                #theta=0.5
                #half_collapse = theta*tsing.tsing_core[core_id]+(1-theta)*tsing.tend_core[core_id]
                #frame_mask[get_time_index(half_collapse)]=True
                self.titles=[ r'$t=0$', r'$t=0.25 t_{\rm{sing}}$', r'$t=0.5 t_{\rm{sing}}$', r'$t=0.75 t_{\rm{sing}}$',\
                    r'$t=t_{\rm{sing}}$', r'$t=t_{\rm{song}}$']
            
            frame_list=thtr.frames[frame_mask]
            rm = rainbow_map(len(frame_list))
            if do_plots:
                axes[0][0].set_title(r'$t=0$')
                axes[0][1].set_title(r'$t=0.25 t_{\rm{sing}}$')
                axes[0][2].set_title(r'$t=0.5 t_{\rm{sing}}$')
                axes[0][3].set_title(r'$t=0.75 t_{\rm{sing}}$')
                axes[0][4].set_title(r'$t=t_{\rm{sing}}$')
                axes[0][5].set_title(r'$t=t_{\rm{end}}$')

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


                args = {'linewidth':0.2, 'c':[0.5]*4}
                ext[-1](RR_sort)
                if do_plots:
                    axes[0][nframe].plot( RR_sort, colors.density_units* M_cuml/V_cuml, **args)
                    ext[0](colors.density_units*M_cuml/V_cuml)
                self.profiles[core_id][frame]['rho']=colors.density_units* M_cuml/V_cuml
                self.profiles[core_id][frame]['V_cuml']=V_cuml


                vr = scrub.vr_rel
                vr_cumsum = np.cumsum( vr[ORDER]*dv_sort)/V_cuml
                #vr_cumsum = np.abs(np.cumsum( vr[ORDER]*dv_sort)/V_cuml)
                #vr_cumsum = np.abs(np.cumsum( rho_sort*vr[ORDER]*dv_sort)/V_cuml)
                self.profiles[core_id][frame]['R']=RR_sort
                self.profiles[core_id][frame]['vr_cumsum']=vr_cumsum
                if do_plots:
                    axes[1][nframe].plot(RR_sort, vr_cumsum, **args)
                    ext[1](vr_cumsum)

                vt = scrub.vt_rel
                vt_cumsum = np.cumsum( vt[ORDER]*dv_sort)/V_cuml
                self.profiles[core_id][frame]['vt_cumsum']=vt_cumsum
                if do_plots:
                    axes[2][nframe].plot(RR_sort, vt_cumsum, **args)
                    ext[2](vt_cumsum)

                if do_plots:
                    axes[3][nframe].plot(RR_sort,EK_cuml/EG_cuml,**args)
                    ext[3](EK_cuml/EG_cuml)
                self.profiles[core_id][frame]['energy']=EK_cuml/EG_cuml
                self.profiles[core_id][frame]['EGmean']=EG_cuml
                self.profiles[core_id][frame]['EKmean']=EK_cuml


                flux = -( DD[ORDER]*vr[ORDER]*dv[ORDER])
                r_bins = np.geomspace( 2e-4, 32/128, 32)
                rcen = 0.5*(r_bins[1:]+r_bins[:-1])
                digitized = np.digitize( RR_sort, r_bins)
                mean_flux  =nar([ flux[ digitized == i].sum() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])
                mass_quant  =nar([ M_cuml[ digitized == i].mean() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])
                #M_cuml = np.cumsum( DD[ORDER]*dv[ORDER])
                ok = ~np.isnan(mean_flux)
                self.profiles[core_id][frame]['rbins']=rcen[ok]
                self.profiles[core_id][frame]['mdot']=mean_flux[ok]/mass_quant[ok]*colors.tff
                self.profiles[core_id][frame]['mass_quant']=mass_quant[ok]

        if do_plots:
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
    Nplots=4
    Ntimes=len(obj.titles)
    fig,axes = plt.subplots(Nplots,Ntimes, figsize=(12,12))
    fig.subplots_adjust(hspace=0,wspace=0)
    for title,ax in zip(obj.titles, axes[0]):
        ax.set_title(title)
    row_dict={'EGmean':0,'EKmean':1,'energy':2, 'mdot':3}
    
    ext = [extents() for n in range(Nplots+1)]
    args = {'linewidth':0.2, 'c':[0.5]*4}
    profs=['EGmean','EKmean','energy','mdot']
    #profs=['energy']
    for nprofile, profile in enumerate(profs):
        print(profile)
        core_list=list(obj.profiles.keys())
        #core_list=core_list[:5]
        for core_id in core_list:
            row = row_dict[profile]
            frames = sorted(list(obj.profiles[core_id].keys()))
            for nframe,frame in enumerate(frames):
                ax = axes[row][nframe]
                if profile != 'mdot':
                    R = nar(obj.profiles[core_id][frame]['R'].v)
                    y = nar(obj.profiles[core_id][frame][profile].v)
                else:
                    R = nar(obj.profiles[core_id][frame]['rbins'])
                    y = nar(obj.profiles[core_id][frame]['mdot'])
                if profile == 'EKmean' and False:
                    #smooth
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
    for ax in axes[3]:
        ax.set(xscale='log',yscale='linear',ylim=ext[3].minmax, ylabel=r'$Dot Log M$', xlim=ext[-1].minmax)
    for row in axes:
        for ax in row[1:]:
            ax.set(ylabel='', yticks=[])
    for ax in axes[-1]:
        ax.set(xlabel='r')
    print('savefig')
    fig.savefig('plots_to_sort/energies.png')

def density(obj):
    Nframes=len(obj.titles)
    Nrows=3
    fig,axes=plt.subplots(Nrows,Nframes)
    fig.subplots_adjust(hspace=0,wspace=0)
    core_list=list(obj.profiles.keys())
    #core_list=core_list[:2] #kludge
    ext=[extents(),extents(),extents()]
    args = {'linewidth':0.2, 'c':[0.5]*4}
    profile='rho'
    #profile='V_cuml'
    #profile = 'mass_quant'
    def nfw(r, rho0, rs):
        return rho0/(r/rs*(1+r/rs)**2)
    profile_list=['rho','mass','mass_over_r']
    for nprof, profile in enumerate(profile_list):
        for core_id in core_list:
            print('fitting',core_id,nprof)
            frames = sorted(list(obj.profiles[core_id].keys()))
            for nframe,frame in enumerate(frames):
                ax=axes[nprof][nframe]
                if profile=='rho':
                    R = nar(obj.profiles[core_id][frame]['R'].v)
                    y = nar(obj.profiles[core_id][frame]['rho'].v) #*R**(-3)
                    V = nar(obj.profiles[core_id][frame]['V_cuml'].v) #*R**(-3)
                    m = y*V
                    y=m*R**(-3)#*R**(-3)
                    ext[0](y)
                if profile=='mass':
                    R = nar(obj.profiles[core_id][frame]['R'].v)
                    y = nar(obj.profiles[core_id][frame]['rho'].v) #*R**(-3)
                    V = nar(obj.profiles[core_id][frame]['V_cuml'].v) #*R**(-3)
                    m = y*V
                    y=m#*R**(-3)
                    ext[1](y)
                if profile=='mass_over_r':
                    R = nar(obj.profiles[core_id][frame]['R'].v)
                    y = nar(obj.profiles[core_id][frame]['rho'].v) #*R**(-3)
                    V = nar(obj.profiles[core_id][frame]['V_cuml'].v) #*R**(-3)
                    m = y*V
                    y=m/R#*R**(-3)
                    ext[2](y)

                if 0:
                    R = nar(obj.profiles[core_id][frame]['R'].v)
                    y = nar(obj.profiles[core_id][frame]['rho'].v) #*R**(-3)
                    V = nar(obj.profiles[core_id][frame]['V_cuml'].v) #*R**(-3)
                    m = y*V
                    y=m*R**(-3)#*R**(-3)
                    #y = V
                if 0:
                    R = nar(obj.profiles[core_id][frame]['R'].v)
                    y = nar(obj.profiles[core_id][frame][profile].v) #*R**(-3)

                if 0:
                    R = nar(obj.profiles[core_id][frame]['rbins'])
                    y = nar(obj.profiles[core_id][frame][profile])

                if 0:
                    try:
                        popt,pcov=curve_fit( nfw, R, y)
                        ax.plot(R, nfw(R,*popt))
                    except:
                        print('bad')

                ext[-1](R)
                ax.plot( R, y, **args)
    for nr,row in enumerate(axes):
        for ax in row:
            ax.set(yscale='log',xscale='log', ylim=ext[nr].minmax, xlim=ext[-1].minmax)
    for row in axes:
        for ax in row[1:]:
            ax.set(ylabel='', yticks=[])
    for ax in axes[-1]:
        ax.set(xlabel='r')
    fig.savefig('plots_to_sort/density_%s'%(obj.this_looper.sim_name))
            
def replotter(obj):
    Nplots=5
    Ntimes=len(obj.titles)
    fig,axes = plt.subplots(Nplots,Ntimes, figsize=(12,12))
    fig.subplots_adjust(hspace=0,wspace=0)
    for title,ax in zip(obj.titles, axes[0]):
        ax.set_title(title)

    
    ext = [extents() for n in range(Nplots+1)]

    args = {'linewidth':0.2, 'c':[0.5]*4}
    #profs=['rho','M/r','vr_cumsum','vt_cumsum','energy']
    profs=['rho','vr_cumsum','vt_cumsum','energy','M/r']
    row_dict={'rho':0,'vr_cumsum':2,'vt_cumsum':3,'energy':4, 'M/r':1}
    #profs=['rho']#, 'energy']
    fit_a=[]
    fit_p=[]
    for nprofile, profile in enumerate(profs):
        print(profile)

        core_list=list(obj.profiles.keys())
        #core_list=core_list[:5] #kludge
        for core_id in core_list:
            row = row_dict[profile]
            frames = sorted(list(obj.profiles[core_id].keys()))
            for nframe,frame in enumerate(frames):
                ax = axes[row][nframe]
                AllR = nar(obj.profiles[core_id][frame]['R'].v)
                if profile == 'M/r':
                    Q = (nar(obj.profiles[core_id][frame]['rho'].v)*\
                        nar(obj.profiles[core_id][frame]['V_cuml'])/\
                        nar(obj.profiles[core_id][frame]['R']))
                else:
                    Q = nar(obj.profiles[core_id][frame][profile].v)
                rbins = nar(obj.profiles[core_id][frame]['rbins'])
                rcen = 0.5*(rbins[1:]+rbins[:-1])
                digitized = np.digitize( AllR, rbins)
                quan = nar([Q[digitized==i].mean() if (digitized==i).any() else np.nan for i in range(1,len(rbins))])
                ok = ~np.isnan(quan)
                y = quan[ok]
                R = rcen[ok]
                ext[-1](R)


                if profile=='rho' and False:
                    if 0:
                        R = nar(obj.profiles[core_id][frame]['R'].v)
                        y = nar(obj.profiles[core_id][frame]['rho'].v) #*R**(-3)
                        V = nar(obj.profiles[core_id][frame]['V_cuml'].v) #*R**(-3)
                        m = y*V
                        y=m*R**(-3)#*R**(-3)
                    if 1:
                        R = nar(obj.profiles[core_id][frame]['rbins'])
                        m = nar(obj.profiles[core_id][frame]['mass_quant'])
                        y = m/R**3

                if profile=='energy':
                    #def fun(R,A,P):
                    #    return nar(A/(R/1e-3)**2)+P
                    y=1/y

                ext[row](y)
                ax.plot( R, y, **args)
                if profile=='energy':

                    if nframe==len(frames)-1:
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
    row = row_dict['rho']
    for nax,ax in enumerate(axes[row]):
        ax.set(xscale='log',yscale='log',ylim=ext[row].minmax, ylabel=r'$<\rho>(<r) [\rm{g/cc}]$', xlim=ext[-1].minmax)
        if nax == len(frames)-1:
            r = np.geomspace(1e-2,0.2)
            ax.plot( r,1e6*(r/0.1)**-2,c='r')

    if 1:
        row = row_dict['M/r']
        for ax in axes[row]:
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax, ylabel=r'$M(<r)/r$',xlim=ext[-1].minmax)
        
    if 1:
        row = row_dict['vr_cumsum']
        for ax in axes[row]:
            ax.set(xscale='log',yscale='linear',ylim=ext[row].minmax, ylabel=r'$<v_r>(<r) [c_s]$', xlim=ext[-1].minmax)
            ax.axhline(0,c=[0.5]*4)

        row = row_dict['vt_cumsum']
        for nframe,ax in enumerate(axes[row]):
            ax.set(xscale='log',yscale='linear',ylim=ext[row].minmax, ylabel=r'$<v_t>(<r)$ [c_s]', xlim=ext[-1].minmax)
            if  nframe == len(frames)-2:
                RR_sort = np.geomspace(1e-4,1e1)
                #ax.plot(RR_sort, RR_sort*vt_cumsum[-1]/RR_sort[-1], c='r')
                ax.plot(RR_sort, 4*(RR_sort/0.1)**0.5, c='r')

        row = row_dict['energy']
        for nframe,ax in enumerate(axes[row]):
            #THIS ONE
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax, ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)

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

            if nframe==len(frames)-1:
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


    for row in axes:
        for ax in row[1:]:
            ax.set(ylabel='', yticks=[])
    for ax in axes[-1]:
        ax.set(xlabel='r')
    print('saving')
    timescale = ['0_tsing','tsing_tsong'][obj.timescale]
    fig.savefig('plots_to_sort/radial_profile_horse_%s_%s'%(obj.this_looper.sim_name,timescale))




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
        timescale = 0 #0= 0-tsing, 1=tsing-tsing
        mp.run(core_list=core_list,tsing=tsing_tool[sim], timescale=timescale)#, r_inflection=anne.inflection[sim])
if 0:
    density(mp)
if 1:
    replotter(mp)
if 0:
    geplotter(mp)
