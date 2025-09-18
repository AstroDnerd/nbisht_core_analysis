

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
                    y = y #* colors.density_units
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

