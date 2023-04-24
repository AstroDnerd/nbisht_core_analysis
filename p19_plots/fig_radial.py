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
#sim_list=['u502']
plt.close('all')
if 0:
    for sim in sim_list:
        loop = TL.loops[sim]
        camera = camera_path.camera_1( loop, 'smooth_zoom_2')
        core_proj_three.core_proj_multiple(loop,axis_list=[0],core_list=[74],frame_list=[0,10],camera=camera, main_core=74)


def gtoy(obj,suffix=''):
    #Nplots=5
    if 0:
        #extra plots for learning
        profs=['rho','vr_cumsum','vt_cumsum','energy','M/r']
        row_dict={'rho':0,'vr_cumsum':2,'vt_cumsum':3,'energy':4, 'M/r':1}
    if 1:
        #the paper version
        row_dict={'rho':0,'vr_cumsum':1,'vt_cumsum':2,'energy':3}
        profs=['rho','vr_cumsum','vt_cumsum','energy']
    if 0:
        #kludge for dev
        row_dict={'rho':0,'vr_cumsum':1,'vt_cumsum':2,'energy':3}
        profs=['vr_cumsum','vt_cumsum']#,'vt_cumsum','energy'
    if 0:
        #kludge for dev
        row_dict={'rho':0,'vr_cumsum':1,'vt_cumsum':2,'energy':0}
        profs=['energy']
    Nplots=len(profs)
    Nplots=4
    Ntimes=len(obj.titles)
    #fig,axes = plt.subplots(Nplots,Ntimes, figsize=(12,12))
    #fig.subplots_adjust(hspace=0,wspace=0)
    fig,ax=plt.subplots(1,1,figsize=(12,12))
    #for title,ax in zip(obj.titles, axes[0]):
    #    ax.set_title(title)

    
    ext = [extents() for n in range(Nplots+1)]

    args = {'linewidth':0.2, 'c':[0.5]*3}

    #profs=['rho']#, 'energy']
    fit_a=[]
    fit_p=[]
    fit_x=[]
    fit_vr=[]
    fit_y_r_eng=[]
    for nprofile, profile in enumerate(profs):
        print(profile)

        profiles = obj.profiles_gas
        core_list=list(profiles.keys())
        n = 5
        #core_list=core_list[n:n+5] #kludge
        #core_list=core_list[0:10] #kludge
        temp_core_list=core_list#[50:60]
        bad_cores=[70,71, 109,93, 177,180,183,198, 248,238, 292]
        thing1 = set(temp_core_list)
        thing2 = set(bad_cores)
        print("bad",bad_cores)
        print("len",len(core_list))
        core_list = list( thing1-thing2)
        #core_list=bad_cores
        print(core_list)
        #core_list=core_list[n:n+1] #kludge
        #print('kludge core list')
        rm = rainbow_map(len(core_list))
        for ncore,core_id in enumerate(core_list):
            c = rm(ncore)
            row = row_dict[profile]
            frames = sorted(list(profiles[core_id].keys()))
            for nframe,frame in enumerate(frames):
                #ax = axes[row][nframe]
                AllR = nar(profiles[core_id][frame]['R'].v)
                if profile == 'M/r':
                    Q = (nar(profiles[core_id][frame]['rho'].v)*\
                        nar(profiles[core_id][frame]['V_cuml'])/\
                        nar(profiles[core_id][frame]['R']))
                else:
                    Q = nar(profiles[core_id][frame][profile].v)
                rbins = nar(profiles[core_id][frame]['rbins'])
                rcen = 0.5*(rbins[1:]+rbins[:-1])
                digitized = np.digitize( AllR, rbins)
                quan = nar([Q[digitized==i].mean() if (digitized==i).any() else np.nan for i in range(1,len(rbins))])
                ok = ~np.isnan(quan)
                y = quan[ok]
                R = rcen[ok]
                #ext[-1](R)
                ext[row](nar([0,3]))

                if profile == 'rho':
                    y = y * colors.density_units
                if profile=='energy':
                    #def fun(R,A,P):
                    #    return nar(A/(R/1e-3)**2)+P
                    y=1/y
                    #if nframe==0:
                    #    y=y*R

                if nframe == 3:
                    ax.plot(R,y,**args)
                    ext[row](y)
                    ext[-1](R)
                if nframe == 3 and False:
                    mmm = y.max()
                    ind = np.argmax(y)
                    if mmm > 10:
                        continue
                    dy = y[1:]-y[:-1]
                    yc = 0.5*(y[1:]+y[:-1])
                    rc = 0.5*(R[1:]+R[:-1])
                    off = 5
                    if (dy<=0).any():
                        first_max = np.where(dy[off:]<=0)[0][0]
                    else:
                        first_max=-1
                    first_max+=off-1
                    mmm=yc.max()
                    ind = np.argmax(yc)
                    Rmax = rc[first_max]
                    y_first=yc[first_max]

                    #plt.scatter( Rmax, yc.max(), marker='*')
                    ext[row](y)
                    RforMax=1
                    R2=rc*RforMax/Rmax
                    y2=yc/y_first
                    #R2=R
                    #ax.plot( R2, y, **args)

                    Q = 0.7#*y_first
                    RQ = R2[ np.argmin( np.abs(y2[:first_max]-Q))]
                    if 0:
                        ax.plot( R2, y2)
                        ax.axvline(RforMax)
                        ext[-1](R2)
                        ax.scatter(RQ,Q)
                        ext[-1](R2)
                    if 1:
                        R3 = R2**(RQ/0.7)
                        ax.plot( R3, y2, c=c)
                        ax.axvline(RforMax)
                        ext[-1](R3)

                    if 1:
                        plt.text( R3[-1],y2[-1],"%s"%core_id)
                    #ext[-1](RQ)
                    print(R2)
                    ext[row](y2)
                    #ext[-1](nar([0.1,100]))


                    #axes[row-1][nframe].plot( rc,dy)
                    #axes[row-1][nframe].axhline(0)




    if 'energy' in profs:
        row = row_dict['energy']
        row=0
        #ax=axes[2][3]
        ax.set(xscale='log',yscale='linear', ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)
        axes=[[ax]]
        for nframe,ax in enumerate(axes[row]):
            #THIS ONE
            #ax.set(xscale='log',yscale='log',ylim=ext[row].minmax, ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)
            ax.set(xscale='log',yscale='linear',ylim=[0,3], ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)

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

            if nframe==0:
                r2 = np.geomspace(1e-2,0.1,32)
                y = 10*(r2/0.1)**-1
                ax.plot(r2,y,c='r')

            if nframe==len(frames)-1 and False:
                RR_sort = np.geomspace(1e-4,1e1)
                mean_a = nar(fit_a).mean()
                mean_p = nar(fit_p).mean()
                mean_a = 0.5
                mean_p=-0.5
                ax.plot(rbins, mean_a*np.log(rbins/1e-3)+mean_p, c='r')

            
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
    timescale = ['0_tsing','tsing_tsung','16_panel'][obj.timescale]
    fig.savefig('plots_to_sort/gtoy_%s_%s_%s'%(obj.this_looper.sim_name,timescale,suffix))

            
def replotter(obj,suffix=''):
    #Nplots=5
    if 0:
        #extra plots for learning
        profs=['rho','vr_cumsum','vt_cumsum','energy','M/r']
        row_dict={'rho':0,'vr_cumsum':2,'vt_cumsum':3,'energy':4, 'M/r':1}
    if 1:
        #the paper version
        row_dict={'rho':0,'vr_cumsum':1,'vt_cumsum':2,'energy':3}
        profs=['rho','vr_cumsum','vt_cumsum','energy']
    if 0:
        #kludge for dev
        row_dict={'rho':0,'vr_cumsum':1,'vt_cumsum':2,'energy':3}
        profs=['vr_cumsum','vt_cumsum']#,'vt_cumsum','energy'
    if 0:
        #kludge for dev
        row_dict={'rho':0,'vr_cumsum':1,'vt_cumsum':2,'energy':3}
        profs=['energy']
    Nplots=len(profs)
    Nplots=4
    Ntimes=len(obj.titles)
    fig,axes = plt.subplots(Nplots,Ntimes, figsize=(12,12))
    fig.subplots_adjust(hspace=0,wspace=0)
    for title,ax in zip(obj.titles, axes[0]):
        ax.set_title(title)

    
    ext = [extents() for n in range(Nplots+1)]

    args = {'linewidth':0.2, 'c':[0.5]*3}

    #profs=['rho']#, 'energy']
    fit_a=[]
    fit_p=[]
    fit_x=[]
    fit_vr=[]
    fit_y_r_eng=[]
    for nprofile, profile in enumerate(profs):
        print(profile)

        profiles = obj.profiles_gas
        core_list=list(profiles.keys())
        #core_list=core_list[:5] #kludge
        #print('kludge core list')
        for core_id in core_list:
            row = row_dict[profile]
            frames = sorted(list(profiles[core_id].keys()))
            for nframe,frame in enumerate(frames):
                ax = axes[row][nframe]
                AllR = nar(profiles[core_id][frame]['R'].v)
                if profile == 'M/r':
                    Q = (nar(profiles[core_id][frame]['rho'].v)*\
                        nar(profiles[core_id][frame]['V_cuml'])/\
                        nar(profiles[core_id][frame]['R']))
                else:
                    Q = nar(profiles[core_id][frame][profile].v)
                rbins = nar(profiles[core_id][frame]['rbins'])
                rcen = 0.5*(rbins[1:]+rbins[:-1])
                digitized = np.digitize( AllR, rbins)
                quan = nar([Q[digitized==i].mean() if (digitized==i).any() else np.nan for i in range(1,len(rbins))])
                ok = ~np.isnan(quan)
                y = quan[ok]
                R = rcen[ok]
                ext[-1](R)

                if profile == 'rho':
                    y = y * colors.density_units
                if profile=='energy':
                    #def fun(R,A,P):
                    #    return nar(A/(R/1e-3)**2)+P
                    ok = y>0
                    y = y[ok]
                    R = R[ok]
                    y=1/y
                    if np.isnan(y).any() or np.isinf(y).any():
                        pdb.set_trace()
                    #if nframe==0:
                    #    y=y*R
                    if nframe==3 and False:
                        fidr=1e-3
                        ind = np.argmin( np.abs(R-fidr))
                        if ind==0:
                            pdb.set_trace()
                        y = y * 0.5/y[ind]

                        #fid2=2e-2
                        #ind2=np.argmin(np.abs(R-fid2))
                        ind2 = np.argmin( np.abs( y-2.01))
                        print("%5.1f %5.1f %5.1f %5.1f"%(ind,ind2, y[ind],y[ind2]))
                        fit_y_r_eng.append(y[ind2])


                ext[row](y)
                ax.plot( R, y, **args)
                if profile == 'vr_cumsum':
                    if nframe==2:
                        pfit = np.polyfit( np.log(R), y,1)
                        fit_vr.append(pfit)
                        #print(fit_vr)
                        #ax.plot( R, pfit[0]*np.log(R)+pfit[1])

                if profile=='energy':
                    if nframe==0:
                        pfit = np.polyfit( np.log(R), np.log(y),1)
                        #print(pfit)
                        fit_x.append(pfit[0])

                    if nframe==len(frames)-1:
                        rmin=1e-3
                        def fun(R,A,P):
                            return A*np.log(R/rmin)+P
                        ok = (R>1e-3)*(R<1e-2)
                        rbins = np.geomspace(5e-4,1e-2,32)
                        if y[ok][0]<2:
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

    if 0:
        def cumr(arr):
            x = nar(arr)+0
            x.sort()
            y=np.arange(x.size)/x.size
            return x,y
        fig2,ax2=plt.subplots(2,2)
        a0=ax2[0][0];a1=ax2[0][1]; a2=ax2[1][0];ax3=ax2[1][1]
        a0.hist(fit_a)
        a1.hist(fit_p)
        a2.hist(fit_x)
        vr_slope = [p[0] for p in fit_vr]
        vr_off  = [p[1] for p in fit_vr]
        xs,ys=cumr(vr_slope)
        ax3.plot(xs,ys,c='r',label='vr slope')
        xs,ys=cumr(vr_off)
        ax3.plot(xs,ys,c='g',label='vr off')


        mean_slope=nar(fit_x).mean()
        a2.axvline(mean_slope)
        a2.axvline(-1.25)
        print(mean_slope)
        fig2.savefig('plots_to_sort/a_p_x_%s.png'%obj.this_looper.sim_name)



    #print(fit_y_r_eng)
    #print(nar(fit_y_r_eng).mean())
    print('set')
    row = row_dict['rho']
    if 'rho' in profs:
        for nax,ax in enumerate(axes[row]):
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax, ylabel=r'$<\rho>(<r) [\rm{g/cc}]$', xlim=ext[-1].minmax)
            if nax == len(frames)-1:
                r = np.geomspace(1e-2,0.2)
                this_y=1e6*(r/0.1)**-2
                ax.plot( r,this_y,c='r')

    if 'M/r' in profs:
        row = row_dict['M/r']
        for ax in axes[row]:
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax, ylabel=r'$M(<r)/r$',xlim=ext[-1].minmax)
        
    if 'vr_cumsum' in profs:
        row = row_dict['vr_cumsum']
        fit_vr = nar(fit_vr)
        #mean_r0=fit_vr[:,0].mean()
        #mean_v0=fit_vr[:,1].mean()
        #print('ug',fit_vr.shape, fit_vr[:,0].shape)
        #print('ug',mean_r0,mean_v0)
        for nax,ax in enumerate(axes[row]):
            ax.set(xscale='log',yscale='linear',ylim=ext[row].minmax, ylabel=r'$<v_r>(<r) [c_s]$', xlim=ext[-1].minmax)
            ax.axhline(0,c=[0.5]*4)
            rbins=np.geomspace(1e-3,1e-2)
            if nax == 2:
                #imagine=-4-np.log10(rbins)
                fitted = -0.5*np.log(rbins/1e-3)-1
                ax.plot(rbins,fitted,c='r')
                #ax.plot(rbins,imagine,c='r')
            if nax==3:
                imagine=np.log(rbins/1e-3)-1
                ax.plot(rbins,imagine,c='r')
                #print(rrr)



    if 'vt_cumsum' in profs:
        row = row_dict['vt_cumsum']
        for nframe,ax in enumerate(axes[row]):
            ax.set(xscale='log',yscale='linear',ylim=ext[row].minmax, ylabel=r'$<v_t>(<r)$ [c_s]', xlim=ext[-1].minmax)
            if  nframe == len(frames)-2 or True:
                RR_sort = np.geomspace(1e-4,1e1)
                #ax.plot(RR_sort, RR_sort*vt_cumsum[-1]/RR_sort[-1], c='r')
                ax.plot(RR_sort, 3*(RR_sort/0.1)**0.5, c='r')

    if 'energy' in profs:
        row = row_dict['energy']
        for nframe,ax in enumerate(axes[row]):
            #THIS ONE
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax, ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)
            #ax.set(xscale='log',yscale='linear',ylim=[0,3], ylabel=r'$EG(<r)/EK(<r)$', xlim=ext[-1].minmax)

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

            if nframe==0:
                r2 = np.geomspace(1e-2,0.1,32)
                y = 10*(r2/0.1)**-1
                ax.plot(r2,y,c='r')

            if nframe==len(frames)-1 and False:
                RR_sort = np.geomspace(1e-4,1e1)
                mean_a = nar(fit_a).mean()
                mean_p = nar(fit_p).mean()
                mean_a = 0.5
                mean_p=-0.5
                ax.plot(rbins, mean_a*np.log(rbins/1e-3)+mean_p, c='r')

            
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
    timescale = ['0_tsing','tsing_tsung','16_panel'][obj.timescale]
    fig.savefig('plots_to_sort/radial_profile_%s_%s_%s'%(obj.this_looper.sim_name,timescale,suffix))




import tsing
reload(tsing)
sim_list=['u501','u502','u503']
if 'tsing_tool' not in dir():
    tsing_tool={}
    for ns,sim in enumerate(sim_list):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()

#import anne
#reload(anne)
##anne.make_inflection()
import radial
reload(radial)
#sim_list0i=['u502']
#mode_list=['Alone']
if 1:
    if 'stuff' not in dir():
        stuff = {}
    sim_list=['u501','u502','u503']
    #sim_list=['u503']
    mode_list=['Alone','Binary','Cluster']
    for sim in sim_list:
        if sim not in stuff:
            stuff[sim]={}
        for mode in mode_list:
            if mode not in stuff[sim]:
                core_list = TL.loops[sim].core_by_mode[mode]
                core_list=core_list[:2]
                thismp=radial.multipro(TL.loops[sim])
                timescale = 2 #0= 0-tsing, 1=tsing-tsing 2=4 panel
                thismp.run(core_list=core_list,tsing=tsing_tool[sim], timescale=timescale,get_particles=False )#, r_inflection=anne.inflection[sim])
                stuff[sim][mode]=thismp
            replotter(stuff[sim][mode],suffix=mode)
if 0:
    sim_list=['u503']
    if 'mp' not in dir():
        for sim in sim_list:
            all_cores=np.unique( TL.loops[sim].tr.core_ids)
            core_list=list(all_cores)
            core_list=None#[323]
            core_list=[323]
            core_list=[25]
            core_list=[74]
            core_list = TL.loops[sim].core_by_mode['Alone']
            core_list=core_list[10:]
            #core_list=None

            mp=radial.multipro(TL.loops[sim])
            timescale = 2 #0= 0-tsing, 1=tsing-tsing 2=4 panel
            mp.run(core_list=core_list,tsing=tsing_tool[sim], timescale=timescale,get_particles=False )#, r_inflection=anne.inflection[sim])
            replotter(mp)
    if 0:
        density(mp)
    if 1:
        replotter(mp)
    if 0:
        geplotter(mp)
    if 0:
        gtoy(mp)
