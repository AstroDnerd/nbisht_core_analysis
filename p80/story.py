
from starter2 import *


from scipy import ndimage

import movie_frames
reload(movie_frames)
import heat_map
reload(heat_map)
plt.close('all')
import pcolormesh_helper as pch
reload(pch)
import progressbar
import time
import eigen
reload(eigen)


import three_loopers_u900 as TL

def plot_monster_6(TTT):
    print('burp')
    core_id=TTT.core_id

    fig,ax=plt.subplots(2,3)
    T0 = (TTT.A0*TTT.Br0*TTT.Bp0)
    T1 = (TTT.A1*TTT.Br1*TTT.Bp1)
    T2 = (TTT.A2*TTT.Br2*TTT.Bp2)
    R = TTT.R

    ax0=ax[0][0]
    ax1=ax[0][1]
    ax2=ax[0][2]
    ax3=ax[1][0]
    ax4=ax[1][1]
    ax5=ax[1][2]

    ext = extents()
    ext( np.abs(T0));ext(np.abs(T1));ext(np.abs(T2))
    big_bins = np.geomspace(ext.minmax[0],ext.minmax[1],64)
    bins = np.linspace(-2,2,128)
    if 1:
        THIS_AX=ax0
        x = np.abs(TTT.R.flatten())
        y = np.abs(TTT.divv.flatten())
        ext=extents()
        ext(x)
        ext(y)
        ext(TTT.rho.flatten())
        bins = np.geomspace(ext.minmax[0],ext.minmax[1],64)

        pch.simple_phase(x,y,ax=THIS_AX,bins=[bins,bins])
        THIS_AX.set(title = 'divv vs R',xlabel='R',ylabel='divv', xscale='log',yscale='log')
    if 1:
        THIS_AX=ax1
        x = np.abs(TTT.rho.flatten())
        y = np.abs(TTT.divv.flatten())
        ext=extents()
        ext(x)
        ext(y)
        #bins = np.geomspace(ext.minmax[0],ext.minmax[1],64)

        pch.simple_phase(x,y,ax=THIS_AX,bins=[bins,bins])
        THIS_AX.set(title = 'divv vs rho', xscale='log',yscale='log')

    if 1:
        THIS_AX=ax2
        x = np.abs(TTT.R.flatten())+0
        x.sort()
        y = np.arange(x.size)/x.size
        THIS_AX.plot(x,y,marker='*')
        THIS_AX.set(xscale='log',yscale='log')
        if 0:
            #this fails
            dy = y[:-1]-y[1:]
            dx = x[:-1]-x[1:]
            yc = 0.5*(y[:-1]+y[1:])
            THIS_AX.plot(yc,dy/dx)
            THIS_AX.set(xscale='log')
    if 1:
        THIS_AX=ax4
        x = np.abs(TTT.R.flatten())+0
        y = np.arange(x.size)/x.size
        x.sort()
        dbins = np.geomspace(x.min(),x.max(),128)
        digitized = np.digitize(x,dbins)
        Cmeans  =nar([ y[ digitized == i].mean() if (digitized==i).any() else np.nan for i in range(0,len(dbins))])

        if 1:
            ok = ~np.isnan(Cmeans)
            x = dbins[ok][1:]
            from scipy.ndimage import gaussian_filter
            y = gaussian_filter(Cmeans[ok][1:],1)

            #THIS_AX.plot(x,y)
            #THIS_AX.set(xscale='log',yscale='log')

        if 1:
            #
            dy = y[:-1]-y[1:]
            dx = x[:-1]-x[1:]
            xc = 0.5*(x[:-1]+x[1:])
            yc = 0.5*(y[:-1]+y[1:])
            THIS_AX.plot(xc,dy/dx)
            #THIS_AX.plot(x,y)
            THIS_AX.set(xscale='log', yscale='log')


        

    if 1:
        THIS_AX=ax3
        x = np.abs(TTT.R.flatten())
        y = np.abs(TTT.rho.flatten())
        ext=extents()
        ext(x)
        ext(y)
        #bins = np.geomspace(ext.minmax[0],ext.minmax[1],64)

        pch.simple_phase(x,y,ax=THIS_AX,bins=[bins,bins])
        THIS_AX.set(title = 'rho vs R', xscale='log',yscale='log')

    fig.savefig('plots_to_sort/wheres_the_beef_c%04d.png'%TTT.core_id)

def plot_monster_7(TTT):
    print('burp')
    core_id=TTT.core_id

    fig,ax=plt.subplots(3,3)

    ext=extents()
    ext(np.abs(TTT.divv))
    big_bins = np.geomspace(ext.minmax[0],ext.minmax[1],64)
    if 1:
        for ni, vi in  enumerate( TTT.dixj):
            for nj, dx in enumerate(vi):
                print(dx.shape)
                THIS_AX=ax[ni][nj]
                divv = np.abs(np.abs(TTT.divv).flatten())
                djvi = np.abs(np.abs(dx).flatten())
                x=djvi
                y=djvi/divv
                #y=djvi/np.abs(TTT.dixj[0][0]).flatten()

                ext=extents()
                ext(x)
                ext(y)
                bins = np.geomspace(ext.minmax[0],ext.minmax[1],64)

                pch.simple_phase(x,y,ax=THIS_AX,bins=[bins,bins])
                THIS_AX.set(xscale='log',yscale='log')

    fig.savefig('plots_to_sort/grads.png')

def meanie(Q,TTT):
    F = np.abs(TTT.R)<5
    Qbar=(TTT.rho*TTT.dv*F*Q).sum()/(TTT.rho*TTT.dv*F).sum()
    return Qbar
def plot_monster_8(objs,sim):

    stuff=defaultdict(list)

    fig,ax=plt.subplots(3,3, figsize=(12,12))
    ax0=ax[0][0]
    ax1=ax[0][1]
    ax2=ax[0][2]
    ax3=ax[1][0]
    ax4=ax[1][1]
    ax5=ax[1][2]
    ax6=ax[2][0]
    ax7=ax[2][1]
    ax8=ax[2][2]
    nintefive=[]
    if 0:
        for core_id in objs:
            TTT = objs[core_id]
            T0 = (TTT.A0*TTT.Br0*TTT.Bp0)
            T1 = (TTT.A1*TTT.Br1*TTT.Bp1)
            T2 = (TTT.A2*TTT.Br2*TTT.Bp2)
            R = TTT.R
            rf = np.abs(R.flatten())+0
            rf.sort()
            cdf = np.arange(rf.size)/rf.size
            ind95=np.where(cdf>0.95)[0][0]
            r95=rf[ind95]
            nintefive.append(r95)

            THIS_AX=ax0
            THIS_AX.plot(rf,cdf)
        THIS_AX.set(yscale='linear',xscale='linear',title='cdf(r)', xlim=[0.1,2.5], ylim=[0.05,1.5])
        THIS_AX.axhline(0.95)
    if 0:
        THIS_AX=ax1
        x = np.array(nintefive)
        x.sort()
        y=np.arange(x.size)/x.size
        THIS_AX.plot(x,y,title='95% of R')

    for core_id in objs:
        TTT=objs[core_id]
        T0 = (TTT.A0*TTT.Br0*TTT.Bp0)
        T1 = (TTT.A1*TTT.Br1*TTT.Bp1)
        T2 = (TTT.A2*TTT.Br2*TTT.Bp2)
        stuff['mean_r'].append(meanie(TTT.R,TTT))
        stuff['mean_T0'].append(meanie(T0,TTT))
        stuff['mean_T1'].append(meanie(T1,TTT))
        stuff['mean_T2'].append(meanie(T2,TTT))
        stuff['mean_A0'].append(meanie(TTT.A0,TTT))
        stuff['mean_A1'].append(meanie(TTT.A1,TTT))
        stuff['mean_A2'].append(meanie(TTT.A2,TTT))
        stuff['mean_Br0'].append(meanie(TTT.Br0,TTT))
        stuff['mean_Br1'].append(meanie(TTT.Br1,TTT))
        stuff['mean_Br2'].append(meanie(TTT.Br2,TTT))
        stuff['mean_Bp0'].append(meanie(TTT.Bp0,TTT))
        stuff['mean_Bp1'].append(meanie(TTT.Bp1,TTT))
        stuff['mean_Bp2'].append(meanie(TTT.Bp2,TTT))
        B0prod = TTT.Br0*TTT.Bp0
        stuff['mean_bprod0'].append(meanie(B0prod,TTT))
        B1prod = TTT.Br1*TTT.Bp1
        stuff['mean_bprod1'].append(meanie(B1prod,TTT))
        B2prod = TTT.Br2*TTT.Bp2
        stuff['mean_bprod2'].append(meanie(B2prod,TTT))
        #stuff['mean_theta'].append(meanie(TTT.theta,TTT))
        also_theta=TTT.Br0*TTT.Bp0+TTT.Br1*TTT.Bp1+TTT.Br2*TTT.Bp2
        stuff['also_theta'].append(meanie(also_theta,TTT))

        stuff['mean_theta_r'].append(meanie(TTT.theta_r,TTT))
        stuff['mean_theta_w'].append(meanie(TTT.theta_w,TTT))
        stuff['mean_theta_b'].append(meanie(TTT.theta_b,TTT))

    for thing in stuff:
        stuff[thing]=nar(stuff[thing])

    if 1:
        THIS_AX=ax0
        x=nar(stuff['mean_r'])
        x.sort()
        y=np.arange(x.size)/x.size
        THIS_AX.plot(x,y)
        THIS_AX.set(title='CDF(R) by core')
    if 1:
        THIS_AX=ax1
        ext=extents()
        x,y=nar(stuff['mean_r']),nar(stuff['mean_T0'])
        ext(x)
        ext(y)
        THIS_AX.scatter(x,y)
        THIS_AX.plot(ext.minmax,ext.minmax)
        THIS_AX.scatter(stuff['mean_r'],stuff['mean_T1'])
        THIS_AX.scatter(stuff['mean_r'],stuff['mean_T2'])
        THIS_AX.set(title='Ti vs R')
    if 1:
        THIS_AX=ax2
        THIS_AX.scatter(nar(stuff['mean_T1'])/nar(stuff['mean_T0']),
                        nar(stuff['mean_T2'])/nar(stuff['mean_T0']))
        THIS_AX.set(title='T2/T0 vs T1/T0')

    if 1:
        THIS_AX=ax3
        x,y=stuff['mean_T0'],stuff['mean_bprod0']
        ext=extents()
        ext(x); ext(y)
        THIS_AX.scatter(x,y)
        x,y=stuff['mean_T1'],stuff['mean_bprod1']
        ext(x); ext(y)
        THIS_AX.scatter(x,y)
        x,y=stuff['mean_T2'],stuff['mean_bprod2']
        ext(x); ext(y)
        THIS_AX.scatter(x,y)
        THIS_AX.plot(ext.minmax,ext.minmax)
        THIS_AX.set(title='<Bi*Ci> vs Ti')
    if 1:
        THIS_AX=ax4
        ext=extents()
        x,y=stuff['mean_T0'],stuff['mean_A0']
        ext(x); ext(y)
        THIS_AX.scatter(x,y)
        x,y=stuff['mean_T1'],stuff['mean_A1']
        ext(x); ext(y)
        THIS_AX.scatter(x,y)
        x,y=stuff['mean_T2'],stuff['mean_A2']
        ext(x); ext(y)
        THIS_AX.scatter(x,y)
        THIS_AX.plot(ext.minmax,ext.minmax)
        THIS_AX.set(title='Ai vs Ti')
    if 1:
        THIS_AX=ax5
        ext=extents()
        #x,y=stuff['mean_T0'],stuff['mean_theta']
        #x,y=stuff['mean_Bp0'],stuff['mean_theta']
        #x,y=stuff['mean_theta'],stuff['also_theta']
        #x,y=stuff['mean_r'],stuff['mean_theta']
        #x,y=stuff['mean_r'],stuff['mean_bprod0']
        #x,y=stuff['mean_r'],stuff['mean_T0']+stuff['mean_T1']+stuff['mean_T2']
        x,y=stuff['mean_r'],stuff['also_theta']
        THIS_AX.scatter(x,y)
        x,y=stuff['mean_r'],stuff['mean_bprod0']
        THIS_AX.scatter(x,y)
        x,y=stuff['mean_r'],stuff['mean_bprod1']
        THIS_AX.scatter(x,y)
        x,y=stuff['mean_r'],stuff['mean_bprod2']
        THIS_AX.scatter(x,y)
        THIS_AX.set(title='BdotC vs R')
        ext(x)
        ext(y)
        THIS_AX.plot(ext.minmax,ext.minmax)

    if 1:
        THIS_AX=ax6
        THIS_AX.scatter(stuff['mean_r'], stuff['mean_theta_r'])
        THIS_AX.set(title='theta R vs R')
    if 1:
        THIS_AX=ax7
        THIS_AX.scatter(stuff['mean_r'], stuff['mean_theta_w'])
        THIS_AX.set(title='Theta W vs R')
    if 1:
        THIS_AX=ax8
        #THIS_AX.scatter(stuff['mean_r'], stuff['mean_theta_b'])
        THIS_AX.scatter(stuff['mean_theta_w'], stuff['mean_theta_b'])
        THIS_AX.set(title='Theta B vs theta W ')

    if 0:
        THIS_AX=ax3
        x,y=nar(stuff['mean_T0']),stuff['mean_Br0']*stuff['mean_Bp0']*stuff['mean_A0']
        #ext(x)
        #ext(y)
        THIS_AX.scatter(x,y)
        #THIS_AX.plot(ext.minmax,ext.minmax)


    fig.tight_layout()
    fig.savefig('plots_to_sort/all_core_%s'%sim)


sim_list=['u902']
new_ddd=False
if 'clobber' not in dir():
    clobber = False
if 'ddd' not in dir() or clobber:
    ddd={}
if 1:
    for sim in sim_list:
        if sim not in ddd:
            ddd[sim]={}

        core_list=None
        #core_list=[7]
        core_list=[74]#, 112]
        #core_list=[112]
        #core_list=[74, 112]
        #core_list=[112]
        core_list=[74]
        core_list=TL.loops[sim].core_by_mode['Alone']
        core_list=core_list

        #core_list=core_list[8:9]
        for core_id in core_list:
            if core_id in ddd[sim]:
                continue
            dddd = eigen.dq_dt2(TL.loops[sim])
            dddd.run(core_list=[core_id], OOM=False)
            ddd[sim][core_id] = dddd

if 1:
    for sim in ddd:
        plot_monster_8(ddd[sim],sim)
if 0:
    for sim in ddd:
        for core_id in ddd[sim]:
            #plot_monster_7(ddd[sim][core_id])
            plot_monster_6(ddd[sim][core_id])
