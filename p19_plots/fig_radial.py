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
import track_loader as TL
sim_list=['u501','u502','u503']
#sim_list=['u502']
plt.close('all')
if 0:
    for sim in sim_list:
        loop = TL.loops[sim]
        camera = camera_path.camera_1( loop, 'smooth_zoom_2')
        core_proj_three.core_proj_multiple(loop,axis_list=[0],core_list=[74],frame_list=[0,10],camera=camera, main_core=74)


            
def replotter(obj,suffix1='', redlines=False, subset=0):
    OneExt = False
    #Nplots=5
    suffix=""
    if subset==0:
        #the paper version
        row_dict={'drho':0,'vr_cumsum':1,'vt_cumsum':2,'energy':3, 'virial':3}
        profs=['drho','vr_cumsum','vt_cumsum','virial']
        #row_dict={'drho':0,'vr_cumsum':1,'vt_cumsum':2,'energy':3, 'virial':3, 'rho_phase':4}
        #profs=['drho','vr_cumsum','vt_cumsum','virial','rho_phase']
        suffix="RVRVTVIR"
    if subset==1:
        #energy plots
        profs=['EGmean', 'EKmean', 'Emag','Etherm', 'virial']
        row_dict={'Emag':2, 'Etherm':3,'EGmean':0, 'EKmean':1, 'virial':4}
        #row_dict={'Emag':0, 'Etherm':0,'EGmean':0, 'EKmean':0}
        suffix="ENG"
        OneExt=True
    if subset==3:
        #energy plots
        profs=['EKmean', 'EGrat', 'ETrat','EBrat']
        row_dict={'EKmean':0, 'EGrat':1, 'ETrat':2, 'EBrat':3}
        #row_dict={'Emag':0, 'Etherm':0,'EGmean':0, 'EKmean':0}
        suffix="ENG"
        OneExt=True
    if subset==2:
        #extra plots for learning
        profs=['rho','EGmean', 'EKmean', 'virial']
        row_dict={'rho':0,'virial':1, 'EGmean':2, 'EKmean':3}
    if subset==4:
        #profs = ['drho', 'dvr', 'dvt', 'virial']
        profs = ['drho', 'dvr', 'dvt', 'virial']
        row_dict = {'drho':0, 'dvr':1,'dvt':2,'virial':3}
    if subset==5:
        #profs = ['drho', 'dvr', 'dvt', 'virial']
        #profs = ['drho', 'pdfrho']
        profs = ['drho', 'rho_phase']
        row_dict = {'drho':0, 'rho_phase':1}

    if subset==2:
        #kludge
        profs=['vr_cumsum']
    if subset==2:
        #field me
        row_dict={'rho':0,'b2':0,'vr_cumsum':1,'vt_cumsum':2,'energy':3}
        profs=['rho','b2']
    if subset==2:
        #kludge for dev
        row_dict={'rho':0,'vr_cumsum':1,'vt_cumsum':2,'energy':3}
        profs=['vr_cumsum','vt_cumsum']#,'vt_cumsum','energy'
    if subset==2:
        #kludge for dev
        row_dict={'rho':0,'vr_cumsum':1,'vt_cumsum':2,'energy':3}
        profs=['energy']
    if subset==2:
        #kludge for dev
        row_dict={'rho':0,'vr_cumsum':1,'britton':2}
        profs=['rho', 'vr_cumsum', 'britton']
    Nplots=len(profs)
    Nplots=max(2, len(profs))
    Ntimes=len(obj.titles)
    fig,axes = plt.subplots(Nplots,Ntimes, figsize=(12,12))
    fig.subplots_adjust(hspace=0,wspace=0)
    for title,ax in zip(obj.titles, axes[0]):
        ax.set_title(title)

    
    if not OneExt:
        ext = [extents() for n in range(Nplots+1)]
    else:
        #this construct also makes a list of extents objects,
        #but uses the same object for the first bunch.
        #ext = [extents()]*Nplots+[extents()]
        ext = [extents()]+3*[extents()]+[extents()]

    args = {'linewidth':0.2, 'c':[0.5]*3}

    #profs=['rho']#, 'energy']
    fit_a=[]
    fit_p=[]
    fit_x=[]
    fit_vr=[]
    fit_y_r_eng=[]
    fig2,axes2=plt.subplots(1,Ntimes, figsize=(12,3))
    extrho=extents()
    for nprofile, profile in enumerate(profs):

        profiles = obj.profiles_gas
        core_list=list(profiles.keys())
        #core_list=core_list[:5] #kludge
        #print('kludge core list')
        for ncore,core_id in enumerate(core_list):
            if core_id == 37:
                print("KLUDGE skipping 37")
                continue
            row = row_dict[profile]
            frames = sorted(list(profiles[core_id].keys()))
            for nframe,frame in enumerate(frames):
                ax = axes[row][nframe]
                AllR = nar(profiles[core_id][frame]['R'].v)
                args['c']=[0.5]*3
                if profile == 'M/r':
                    Q = (nar(profiles[core_id][frame]['rho'].v)*\
                        nar(profiles[core_id][frame]['V_cuml'])/\
                        nar(profiles[core_id][frame]['R']))
                elif profile == 'b2':
                    B = nar(profiles[core_id][frame]['b2'])
                    dv = nar(profiles[core_id][frame]['dv_sort'])
                    Vc = nar(profiles[core_id][frame]['V_cuml'])
                    Q = np.cumsum(B*dv)/Vc
                    args['c']=[1.0,0.0,0.0]
                elif profile == 'britton_tcross':
                    #Q = nar(profiles[core_id][frame]['britton'])
                    sqrtG=np.sqrt(colors.G)
                    G = colors.G
                    rho = nar(profiles[core_id][frame]['rho'])
                    R = nar(profiles[core_id][frame]['R'])
                    R0 = 1e-2

                    ind = np.argmin( np.abs(R-R0))
                    rho0 = rho[ind]
                    R_rho = R0*np.sqrt(rho/rho0)
                    vr  = np.abs(nar(profiles[core_id][frame]['vr_cumsum']))
                    t_cross = R_rho/(vr+1)
                    t_ff = np.sqrt( 3*np.pi/32/(G*rho))

                    Q = t_cross/t_ff
                elif profile == 'britton':
                    M = nar(profiles[core_id][frame]['rho'])*\
                        nar(profiles[core_id][frame]['V_cuml'])
                    Q = M
                    R0=1e-2
                elif profile == 'virial':
                    #Q = nar(profiles[core_id][frame]['energy'])
                    EK=nar(profiles[core_id][frame]['EKmean'])
                    EG=nar(profiles[core_id][frame]['EGmean'])
                    dv = nar(profiles[core_id][frame]['dv_sort'])
                    Vc = nar(profiles[core_id][frame]['V_cuml'])
                    #Q = np.cumsum(B*dv)/Vc
                    EBlocal=nar(profiles[core_id][frame]['b2_sort']**2/np.sqrt(4*np.pi))
                    EB = np.cumsum(EBlocal*dv)/Vc
                    rho=nar(profiles[core_id][frame]['rho_sort'])
                    c=1
                    ETlocal = c**2*rho*np.log(rho)
                    ET =np.abs( np.cumsum(ETlocal*dv)/Vc)
                    #bot += TE
                    #Q = np.abs(np.abs(EK+EB+ET-EG))
                    #Q = EK/EB
                    #Q = 2*EK/(EG+ET+EB)
                    Q = EK/EG
                elif profile == 'Emag':
                    #B = nar(profiles[core_id][frame]['b2_sort'])**(1/0.3)/np.sqrt(4*np.pi)/1000
                    B = nar(profiles[core_id][frame]['b2_sort'])**(2)/np.sqrt(4*np.pi)
                    dv = nar(profiles[core_id][frame]['dv_sort'])
                    Vc = nar(profiles[core_id][frame]['V_cuml'])
                    Q = np.cumsum(B*dv)/Vc
                    args['c']=[0.0,0.0,1.0,1.0]
                elif profile == 'Etherm':
                    rho=nar(profiles[core_id][frame]['rho_sort'])
                    c=1
                    TE = c**2*rho*np.log(rho)
                    dv = nar(profiles[core_id][frame]['dv_sort'])
                    Vc = nar(profiles[core_id][frame]['V_cuml'])
                    Q = np.cumsum(TE*dv)/Vc
                    args['c']=[0.0,1.0,0.0,1.0]
                elif profile == 'EGmean':
                    Q = nar(profiles[core_id][frame][profile].v)
                    args['c'] = [1.0,0.0,0.0,1.0]
                elif profile == 'EKmean':
                    Q = nar(profiles[core_id][frame][profile].v)
                    args['c'] = 'k'# [0.5,0.0,1.0,1.0]
                elif profile == 'EGrat':
                    EG=nar(profiles[core_id][frame]['EGmean'].v)
                    EK = nar(profiles[core_id][frame]['EKmean'].v)
                    Q = EK/EG
                elif profile == 'ETrat':
                    rho=nar(profiles[core_id][frame]['rho_sort'])
                    dv = nar(profiles[core_id][frame]['dv_sort'])
                    Vc = nar(profiles[core_id][frame]['V_cuml'])
                    c=1
                    ETlocal = c**2*rho*np.log(rho)
                    ET =np.abs( np.cumsum(ETlocal*dv)/Vc)
                    EK = nar(profiles[core_id][frame]['EKmean'].v)
                    Q = EK/ET
                elif profile == 'EBrat':
                    B = nar(profiles[core_id][frame]['b2_sort'])**(2)/np.sqrt(4*np.pi)
                    dv = nar(profiles[core_id][frame]['dv_sort'])
                    Vc = nar(profiles[core_id][frame]['V_cuml'])
                    EB = np.cumsum(B*dv)/Vc
                    EK = nar(profiles[core_id][frame]['EKmean'].v)
                    Q = EK/EB
                elif profile == 'pdfrho':
                    import tools.equal_probability_binner as epb
                    reload(epb)
                    rho = profiles[core_id][frame]['drho']+0
                    rho.sort()
                    extrho(rho)
                    y = (np.arange(rho.size)/rho.size)[::-1]
                    ax2=axes2[nframe]
                    ax2.plot(rho,y)
                    #pdf, centers, widths=epb.equal_prob(np.log10(rho.v), 16)
                    #ax2.bar( centers, pdf, width=widths, facecolor=None,edgecolor='k')
                    continue

                elif profile == 'rho_phase':
                    rho = profiles[core_id][frame]['drho']
                    #ok = np.where(rho>1e4)[0]
                    #if len(ok) :
                    #    last_dense = nar(ok).max()

                    #    my_r=AllR[last_dense]*colors.length_units_au
                    #    print(my_r)
                    #    ax.axvline( my_r)
                    #    print('wtf',nframe)
                    ax.axhline( 1e4*colors.density_units)
                    #ax.scatter(AllR*colors.length_units_au,rho*colors.density_units)
                    import pcolormesh_helper as pch
                    reload(pch)
                    X,Y=AllR*colors.length_units_au,rho*colors.density_units
                    pch.simple_phase( X,Y,log=True,ax=ax)

                    ax.set(yscale='log',xscale='log')
                    continue
                else:
                    Q = nar(profiles[core_id][frame][profile])
                    if np.isnan(Q).any():
                        pdb.set_trace()
                rbins = nar(profiles[core_id][frame]['rbins'])
                if len(rbins) < 3:
                    print("NOT ENOUGH BINS")
                    continue
                rcen = 0.5*(rbins[1:]+rbins[:-1])
                digitized = np.digitize( AllR, rbins)
                quan = nar([Q[digitized==i].mean() if (digitized==i).any() else np.nan for i in range(1,len(rbins))])
                ok = ~np.isnan(quan)
                y = quan[ok]
                R = rcen[ok]
                R *= colors.length_units_au

                ext[-1](R)

                if profile == 'rho' or profile == 'drho':
                    y = y * colors.density_units
                    #if (Q>1e4).any():
                    #    dense = Q>1e4
                    #    r_dense = AllR[dense]
                    #    r_pseudo=r_dense.max()*colors.length_units_au
                    #    #ax.axvline(r_pseudo)
                    #    ax.axvline(1200)
                    #    print(r_pseudo)




                ext[row](y)
                ax.plot( R, y, **args)

    #for ax2 in axes2:
    #    ax2.set(yscale='log',xscale='log', xlim=extrho.minmax)
    #fig2.tight_layout()
    #fig2.savefig('plots_to_sort/pdf')
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
        fig2.savefig('plots_to_sort/a_p_x_%s.png'%obj.mon.name)



    #print(fit_y_r_eng)
    #print(nar(fit_y_r_eng).mean())
    for row in axes:
        row[-1].axvline(1200,c=[0.5]*4)
    if 'rho_phase' in profs:
        row = row_dict['rho_phase']
        for nax,ax in enumerate(axes[row]):
            ax.set(xscale='log',yscale='log',ylim=ext[0].minmax, ylabel=r'$\overline{n(r)}~~ [\rm{cm}^{-3}]$', xlim=ext[-1].minmax)
    if 'rho' in profs:
        row = row_dict['rho']
        raise
        for nax,ax in enumerate(axes[row]):
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax, ylabel=r'$\overline{n(<r)}~~ [\rm{cm}^{-3}]$', xlim=ext[-1].minmax)
            if nax == len(frames)-1:
                r = np.geomspace(1e-3,0.2)* colors.length_units_au
                this_y=1e6*(r/1e4)**-(2)
                if redlines:
                    ax.plot( r,this_y,c='r')
    if 'drho' in profs:
        row = row_dict['drho']
        for nax,ax in enumerate(axes[row]):
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax, ylabel=r'$\overline{n(r)}~~ [\rm{cm}^{-3}]$', xlim=ext[-1].minmax)
            if nax == len(frames)-1:
                r = np.geomspace(1e-3,0.2)* colors.length_units_au
                this_y=1e6*(r/1e4)**-(2)
                if redlines:
                    ax.plot( r,this_y,c='r')

    if 'britton' in profs:
        row = row_dict['britton']
        for ax in axes[row]:
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax, ylabel=r'$Tcross/Tff$',xlim=ext[-1].minmax)
            ax.axhline(1)
            ax.axvline(R0)
    if 'b2' in profs:
        row = row_dict['b2']
        for ax in axes[row]:
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax, ylabel=r'$B^2$',xlim=ext[-1].minmax)
    if 'M/r' in profs:
        row = row_dict['M/r']
        for ax in axes[row]:
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax, ylabel=r'$M(<r)/r$',xlim=ext[-1].minmax)
    if 'EGrat' in profs:
        row = row_dict['EGrat']
        for ax in axes[row]:
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax, ylabel=r'$E_K/E_G$',xlim=ext[-1].minmax)
            ax.axhline(0.5, c='k', linewidth=0.1)
            ax.axhline(2, c='k', linewidth=0.1)
    if 'ETrat' in profs:
        row = row_dict['ETrat']
        for ax in axes[row]:
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax, ylabel=r'$E_K/E_T$',xlim=ext[-1].minmax)
            ax.axhline(0.5, c='k', linewidth=0.1)
            ax.axhline(2, c='k', linewidth=0.1)
    if 'EBrat' in profs:
        row = row_dict['EBrat']
        for ax in axes[row]:
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax, ylabel=r'$E_K/E_B$',xlim=ext[-1].minmax)
            ax.axhline(0.5, c='k', linewidth=0.1)
            ax.axhline(2, c='k', linewidth=0.1)
    if 'EGmean' in profs:
        rowG = row_dict['EGmean']
        rowK = row_dict['EKmean']
        ext[rowG](ext[rowK].minmax)
        ext[rowK](ext[rowG].minmax)

    if 'Emag' in profs:
        row = row_dict['Emag']
        for ax in axes[row]:
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax,ylabel='EB(<r>)/V(R)', xlim=ext[-1].minmax)
    if 'Etherm' in profs:
        row = row_dict['Etherm']
        for ax in axes[row]:
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax,ylabel='ET(<r>)/V(R)', xlim=ext[-1].minmax)
    if 'virial' in profs:
        row = row_dict['virial']
        for ax in axes[row]:
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax,ylabel='EK(<r)/EG(<r)', xlim=ext[-1].minmax)
            ax.axhline(0.5,c=[0.5]*4)
            ax.axhline(2.0,c=[0.5]*4)

    if 'EGmean' in profs:
        row = row_dict['EGmean']
        for ax in axes[row]:
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax,ylabel='EG(<r)/V(r)', xlim=ext[-1].minmax)
    if 'EKmean' in profs:
        row = row_dict['EKmean']
        for ax in axes[row]:
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax,ylabel='EK(<r)/V(r)', xlim=ext[-1].minmax)
        
    #if 'dvr' in profs:
    #    row = row_dict['dvr']
    #    fit_vr = nar(fit_vr)
    #    for nax,ax in enumerate(axes[row]):
    #        vr_ylim=ext[row].minmax
    #        vr_ylim = -4,2
    #        ax.set(xscale='log',yscale='linear',ylim=vr_ylim, ylabel=r'$\overline{v_R(<r)}/c_s$', xlim=ext[-1].minmax)
    #if 'dvt' in profs:
    #    row = row_dict['dvt']
    #    for nax,ax in enumerate(axes[row]):
    #        ax.set(xscale='log',yscale='linear',ylim=ext[row].minmax, ylabel=r'$\overline{v_T(<r)}/c_s$', xlim=ext[-1].minmax)
    if 'vr_cumsum' in profs:
        row = row_dict['vr_cumsum']
        fit_vr = nar(fit_vr)
        #mean_r0=fit_vr[:,0].mean()
        #mean_v0=fit_vr[:,1].mean()
        #print('ug',fit_vr.shape, fit_vr[:,0].shape)
        #print('ug',mean_r0,mean_v0)
        for nax,ax in enumerate(axes[row]):
            vr_ylim=ext[row].minmax
            vr_ylim = -4,2
            ax.set(xscale='log',yscale='linear',ylim=vr_ylim, ylabel=r'$\overline{v_R(<r)}/c_s$', xlim=ext[-1].minmax)
            ax.axhline(0,c=[0.5]*4)
            ax.axhline(-1,c=[0.5]*4)
            rbins=np.geomspace(1e-3,1e-2)* colors.length_units_au
            if nax == 2 and redlines:
                #imagine=-4-np.log10(rbins)
                r0=1e4
                fitted = -0.5*(np.log(rbins/r0))
                ax.plot(rbins,fitted,c='r')
                #print(rbins, fitted)
                #ax.plot(rbins,imagine,c='r')
            if nax==3 and redlines:
                imagine=np.log(rbins/1e3)
                ax.plot(rbins,imagine,c='r')
                #print(rrr)



    if 'vt_cumsum' in profs:
        row = row_dict['vt_cumsum']
        for nframe,ax in enumerate(axes[row]):
            ax.set(xscale='log',yscale='linear',ylim=ext[row].minmax, ylabel=r'$\overline{v_T(<r)}/c_s$', xlim=ext[-1].minmax)
            if  nframe == len(frames)-2 and redlines:
                RR_sort = np.geomspace(1e-4,1e1)*colors.length_units_au
                #ax.plot(RR_sort, RR_sort*vt_cumsum[-1]/RR_sort[-1], c='r')
                ax.plot(RR_sort, 9*(RR_sort/colors.length_units_au)**0.5, c='r')

    if 'energy' in profs:
        row = row_dict['energy']
        for nframe,ax in enumerate(axes[row]):
            #THIS ONE
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax, ylabel=r'$\frac{EK(<r)}{EG(<r)}$', xlim=ext[-1].minmax)
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

            if nframe==0 and redlines:
                r2 = np.geomspace(1e-2,0.1,32)
                y = 10*(r2/0.1)**-1
                ax.plot(r2,y,c='r')



    for row in axes:
        for ax in row[1:]:
            ax.set(ylabel='', yticks=[])
    for ax in axes[-1]:
        ax.set(xlabel='R [AU]', xlim=ext[-1].minmax, xscale='log')
    print('saving')
    timescale = ['0_tsing','tsing_tsung','16_panel'][obj.timescale]
    suffix2="%s_%s"%(suffix1,suffix)
    fig.savefig('plots_to_sort/radial_profile_%s_%s_%s.pdf'%(obj.mon.name,timescale,suffix2))
    #fig.savefig('plots_to_sort/radial_profile_%s_%s_%s.png'%(obj.mon.name,timescale,suffix2))




import tsing
reload(tsing)
sim_list=['u501','u502','u503']
TL.load_tracks(sim_list)
if 'tsing_tool' not in dir():
    tsing_tool={}
    for ns,sim in enumerate(sim_list):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()

#import anne
#reload(anne)
##anne.make_inflection()
import radial2
reload(radial2)
#sim_list0i=['u502']
#mode_list=['Alone']
if 1:
    if 'stuff' not in dir():
        stuff = {}
    #sim_list=['u501','u502','u503']
    sim_list=['u501']
    import monster
    reload(monster)
    monster.load(sim_list)
    mode_list=['A']#,'Binary','Cluster']
    for sim in sim_list:
        if sim not in stuff:
            stuff[sim]={}
        for mode in mode_list:
            print("Do ",sim,mode)
            if mode not in stuff[sim]:
                core_list = TL.loops[sim].core_by_mode[mode]
                #core_list=[8]#core_list[:2]
                thismp=radial2.multipro2(monster.closet[sim])
                timescale = 2 #0= 0-tsing, 1=tsing-tsing 2=4 panel
                thismp.run(core_list=core_list, timescale=timescale,get_particles=False, save_sorts=True )#, r_inflection=anne.inflection[sim])
                stuff[sim][mode]=thismp
            replotter(stuff[sim][mode],suffix1=mode, redlines=True, subset=0)

