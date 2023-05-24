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

            
def replotter(obj,suffix='', redlines=False):
    OneExt = False
    #Nplots=5
    if 1:
        #the paper version
        row_dict={'rho':0,'vr_cumsum':1,'vt_cumsum':2,'energy':3, 'virial':3}
        profs=['rho','vr_cumsum','vt_cumsum','virial']
    if 0:
        #energy plots
        profs=['EGmean', 'EKmean', 'Emag','Etherm', 'virial']
        row_dict={'Emag':2, 'Etherm':3,'EGmean':0, 'EKmean':1, 'virial':4}
        #row_dict={'Emag':0, 'Etherm':0,'EGmean':0, 'EKmean':0}
        OneExt=True
    if 0:
        #extra plots for learning
        profs=['rho','EGmean', 'EKmean', 'virial']
        row_dict={'rho':0,'virial':1, 'EGmean':2, 'EKmean':3}

    if 0:
        #kludge
        profs=['vr_cumsum']
    if 0:
        #field me
        row_dict={'rho':0,'b2':0,'vr_cumsum':1,'vt_cumsum':2,'energy':3}
        profs=['rho','b2']
    if 0:
        #kludge for dev
        row_dict={'rho':0,'vr_cumsum':1,'vt_cumsum':2,'energy':3}
        profs=['vr_cumsum','vt_cumsum']#,'vt_cumsum','energy'
    if 0:
        #kludge for dev
        row_dict={'rho':0,'vr_cumsum':1,'vt_cumsum':2,'energy':3}
        profs=['energy']
    if 0:
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
        ext = [extents()]*Nplots+[extents()]

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
                    ET = np.cumsum(ETlocal*dv)/Vc
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
                       


                else:
                    Q = nar(profiles[core_id][frame][profile].v)
                rbins = nar(profiles[core_id][frame]['rbins'])
                rcen = 0.5*(rbins[1:]+rbins[:-1])
                digitized = np.digitize( AllR, rbins)
                quan = nar([Q[digitized==i].mean() if (digitized==i).any() else np.nan for i in range(1,len(rbins))])
                ok = ~np.isnan(quan)
                y = quan[ok]
                R = rcen[ok]
                R *= colors.length_units_au



                if 0:
                    y /= R**(-2.8)
                    #normalize everyone
                    args['c']=[0.5]*4
                    i5000 = np.argmin(np.abs(R-5000))
                    y5000 = y[i5000]
                    y /= 1e3*y5000
                ext[-1](R)

                if profile == 'rho':
                    y = y * colors.density_units
                #if profile=='energy':
                #    ok = y>0
                #    y = y[ok]
                #    R = R[ok]
                #    y=1/y
                #    if np.isnan(y).any() or np.isinf(y).any():
                #        pdb.set_trace()

                ext[row](y)
                ax.plot( R, y, **args)

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
    if 'rho' in profs:
        row = row_dict['rho']
        for nax,ax in enumerate(axes[row]):
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax, ylabel=r'$\overline{\rho(<r)}~~ [\rm{g/cc}]$', xlim=ext[-1].minmax)
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
            ax.axhline(0.5,c=[0.5]*3)
            ax.axhline(2.0,c=[0.5]*3)

    if 'EGmean' in profs:
        row = row_dict['EGmean']
        for ax in axes[row]:
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax,ylabel='EG(<r)/V(r)', xlim=ext[-1].minmax)
    if 'EKmean' in profs:
        row = row_dict['EKmean']
        for ax in axes[row]:
            ax.set(xscale='log',yscale='log',ylim=ext[row].minmax,ylabel='EK(<r)/V(r)', xlim=ext[-1].minmax)
        
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
            ax.set(xscale='log',yscale='linear',ylim=vr_ylim, ylabel=r'$\overline{v_r(<r)}~~[c_s]$', xlim=ext[-1].minmax)
            ax.axhline(0,c=[0.5]*4)
            rbins=np.geomspace(1e-3,1e-2)* colors.length_units_au
            if nax == 2 and redlines:
                #imagine=-4-np.log10(rbins)
                r0=1e4
                fitted = -0.5*(np.log(rbins/r0))
                ax.plot(rbins,fitted,c='r')
                #ax.plot(rbins,imagine,c='r')
            if nax==3 and redlines:
                imagine=np.log(rbins/1e3)
                ax.plot(rbins,imagine,c='r')
                #print(rrr)



    if 'vt_cumsum' in profs:
        row = row_dict['vt_cumsum']
        for nframe,ax in enumerate(axes[row]):
            ax.set(xscale='log',yscale='linear',ylim=ext[row].minmax, ylabel=r'$\overline{v_t(<r)}~~ [c_s]$', xlim=ext[-1].minmax)
            if  nframe == len(frames)-2 and redlines:
                RR_sort = np.geomspace(1e-4,1e1)*colors.length_units_au
                #ax.plot(RR_sort, RR_sort*vt_cumsum[-1]/RR_sort[-1], c='r')
                ax.plot(RR_sort, 3*(RR_sort/1e5)**0.5, c='r')

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
if 0:
    if 'stuff' not in dir():
        stuff = {}
    sim_list=['u501','u502','u503']
    #sim_list=['u503']
    mode_list=['Alone','Binary','Cluster']
    for sim in sim_list:
        if sim not in stuff:
            stuff[sim]={}
        for mode in mode_list:
            print("Do ",sim,mode)
            if mode not in stuff[sim]:
                core_list = TL.loops[sim].core_by_mode[mode]
                #core_list=core_list[:2]
                thismp=radial.multipro(TL.loops[sim])
                timescale = 2 #0= 0-tsing, 1=tsing-tsing 2=4 panel
                thismp.run(core_list=core_list,tsing=tsing_tool[sim], timescale=timescale,get_particles=False )#, r_inflection=anne.inflection[sim])
                stuff[sim][mode]=thismp
            replotter(stuff[sim][mode],suffix=mode)

if 1:
    sim_list=['u501']
    if 'mp' not in dir():
        for sim in sim_list:
            all_cores=np.unique( TL.loops[sim].tr.core_ids)
            core_list=list(all_cores)
            core_list=None#[323]
            core_list=[323]
            core_list=[25]
            core_list=[74]
            core_list = TL.loops[sim].core_by_mode['Alone']
            #core_list=core_list[10:]
            #core_list=None

            mp=radial.multipro(TL.loops[sim])
            timescale = 2 #0= 0-tsing, 1=tsing-tsing 2=4 panel
            mp.run(core_list=core_list,tsing=tsing_tool[sim], timescale=timescale,get_particles=False,save_sorts=True )#, r_inflection=anne.inflection[sim])
            replotter(mp)
    if 1:
        replotter(mp, redlines=True)
    if 0:
        geplotter(mp)
    if 0:
        gtoy(mp)
    if 0:
        density(mp)
