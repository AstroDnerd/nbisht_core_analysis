
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

import tsing
reload(tsing)
import radial
reload(radial)
sim_list=['u502']
if 'tsing_tool' not in dir():
    tsing_tool={}
    for ns,sim in enumerate(sim_list):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()

def phasor_vr_vt(obj):
    Ntimes=len(obj.titles)
    ext = [extents() for n in range(Ntimes+1)]
    row_dict={'rho_hist':0, 'vr_hist':1}
    means_par=[defaultdict(list) for n in range(Ntimes)]
    means_gas=[defaultdict(list) for n in range(Ntimes)]

    profiles = obj.profiles_gas
    core_list=list(profiles.keys())
    #core_list=core_list[:1] #kludge
    for core_id in core_list:
        fig,axes = plt.subplots(4, Ntimes,figsize=(12,12))
        fig.subplots_adjust(wspace=0)
        for title,ax in zip(obj.titles, axes[0]):
            ax.set_title(title)
        frames = sorted(list(profiles[core_id].keys()))
        for nframe,frame in enumerate(frames):
            Pgas = obj.profiles_gas[core_id][frame]
            Ppart= obj.profiles_part[core_id][frame]

            rho_bins=np.geomspace(1e-1,1e8)
            #rho_bins=np.geomspace(1e1,1e5)
            vr_bins=np.linspace(-10,10)
            vt_bins=np.linspace(0,10)


            ax = axes[0][nframe]
            pch.simple_phase(Pgas['vt_sort'], Pgas['vr_sort'], ax=ax, bins=[vt_bins,vr_bins])
            ax.scatter( Ppart['vt_sort'], Ppart['vr_sort'] , s=0.2, c='k')
            ax.set(xscale='linear')

            ax = axes[1][nframe]
            pch.simple_phase(Pgas['rho_sort'], Pgas['vt_sort'], ax=ax, bins=[rho_bins,vt_bins])
            ax.scatter( Ppart['rho_sort'], Ppart['vt_sort'] , s=0.2,c='k')
            ax.set(xscale='log')

            if 0:
                #angular momentum vs rho.
                #its clear from the plot that SPECIFIC angular momentum is the quantity that matters.
                ax = axes[2][nframe]
                #pch.simple_phase(Pgas['rho_sort']*Pgas['vt_sort'],np.abs(Pgas['rho_sort']* Pgas['vr_sort']), ax=ax, bins=[rho_bins,rho_bins])
                #ax.scatter( Ppart['rho_sort']*Ppart['vt_sort'], np.abs(Ppart['rho_sort']*Ppart['vr_sort'] ), s=0.2, c='k')
                pch.simple_phase(Pgas['rho_sort'],np.abs(Pgas['rho_sort']* Pgas['vt_sort']), ax=ax, bins=[rho_bins,rho_bins])
                ax.scatter( Ppart['rho_sort'], np.abs(Ppart['rho_sort']*Ppart['vt_sort'] ), s=0.2, c='k')
                ax.set(xscale='log', yscale='log')

            ax=axes[3][nframe]
            order_gas = np.argsort(Pgas['vt_sort'])
            order_part = np.argsort(Ppart['vt_sort'])
            vtg = Pgas['vt_sort'][order_gas]
            vtp = Ppart['vt_sort'][order_part]

            yg = np.arange(vtg.size)/vtg.size
            yp = np.arange(vtp.size)/vtp.size
            ax.plot( vtg,yg,c='g')
            ax.plot( vtp,yp,c='purple')

            ax=axes[2][nframe]
            order_gas = np.argsort(Pgas['vr_sort'])
            order_part = np.argsort(Ppart['vr_sort'])
            vtg = Pgas['vr_sort'][order_gas]
            vtp = Ppart['vr_sort'][order_part]

            yg = np.arange(vtg.size)/vtg.size
            yp = np.arange(vtp.size)/vtp.size
            ax.plot( vtg,yg,c='g')
            ax.plot( vtp,yp,c='purple')


            #mass_gas = Pgas['rho_sort']*Pgas['dv_sort']
            #mass_part = Ppart['rho_sort']*Ppart['dv_sort']
            #print(mass_gas.sum())
            #mass_cuml_gas = np.cumsum( mass_gas[order_gas])
            #mass_cuml_part = np.cumsum( mass_part[order_part])
            #ax.plot( vtg, mass_cuml_gas)
            #ax.plot( vtp, mass_cuml_part)

        outname='plots_to_sort/phase_frames_%s_c%04d.png'%(obj.this_looper.sim_name, core_id)
        fig.savefig(outname)
        print(outname)

def phasor(obj):
    Ntimes=len(obj.titles)
    ext = [extents() for n in range(Ntimes+1)]
    row_dict={'rho_hist':0, 'vr_hist':1}
    means_par=[defaultdict(list) for n in range(Ntimes)]
    means_gas=[defaultdict(list) for n in range(Ntimes)]

    profiles = obj.profiles_gas
    core_list=list(profiles.keys())
    core_list=core_list[:5] #kludge
    for core_id in core_list:
        fig,axes = plt.subplots(2, Ntimes,figsize=(12,12))
        fig.subplots_adjust(wspace=0)
        for title,ax in zip(obj.titles, axes[0]):
            ax.set_title(title)
        frames = sorted(list(profiles[core_id].keys()))
        for nframe,frame in enumerate(frames):
            Pgas = obj.profiles_gas[core_id][frame]
            Ppart= obj.profiles_part[core_id][frame]

            rho_bins=np.geomspace(1e-1,1e8)
            v_bins=np.linspace(-10,10)


            ax = axes[0][nframe]
            pch.simple_phase(Pgas['rho_sort'], Pgas['vr_sort'], ax=ax, bins=[rho_bins,v_bins])
            ax.set(xscale='log')

            ax = axes[1][nframe]
            pch.simple_phase(Ppart['rho_sort'], Ppart['vr_sort'], ax=ax, bins=[rho_bins,v_bins])
            ax.set(xscale='log')

        fig.savefig('plots_to_sort/phase_frames_%s_c%04d.png'%(obj.this_looper.sim_name, core_id))
if 0:
    if 1:
        row = row_dict['vr_hist']
        bins = np.arange(-8,2,0.5)
        for nframe in range(Ntimes):
            ax=axes[row][nframe]
            #qp = nar(means_par[nframe]['vr'])+0
            #qp.sort()
            #qg = nar(means_gas[nframe]['vr'])+0
            #qg.sort()
            #yp = np.arange( qp.size)/qp.size
            #yg = np.arange( qg.size)/qg.size
            #axes[1][nframe].plot( qp, yp, c='purple')
            #axes[1][nframe].plot( qg, yg, c='g')

            cuml=False
            ax.hist(means_par[nframe]['vr'], histtype='step',color='purple', bins=bins, cumulative=cuml, density=True)
            ax.hist(nar(means_gas[nframe]['vr']).flatten(), histtype='step',color='g', bins=bins, cumulative=cuml, density=True)
    if 'rho_hist' in row_dict:
        row=row_dict['rho_hist']
        for nax, ax in enumerate(axes[row]):
            ax.set(xscale='log',yscale='log',xlabel='rho', ylabel='P(rho)')
    if 'vr_hist' in row_dict:
        row=row_dict['vr_hist']
        for nax, ax in enumerate(axes[row]):
            ax.set(xscale='linear',yscale='log',xlabel='rho', ylabel='P(rho)')
    for row in axes:
        for ax in row[1:]:
            ax.set(ylabel='', yticks=[])

    fig.savefig('plots_to_sort/sub_or_super_%s'%obj.this_looper.sim_name)
    plt.close(fig)


def souper(obj):
    profs = ['rho_hist', 'vr_hist']
    Nplots=len(profs)
    Ntimes=len(obj.titles)
    fig,axes = plt.subplots(Nplots,Ntimes, figsize=(12,12))
    fig.subplots_adjust(wspace=0)
    for title,ax in zip(obj.titles, axes[0]):
        ax.set_title(title)
    ext = [extents() for n in range(Nplots+1)]
    row_dict={'rho_hist':0, 'vr_hist':1}
    means_par=[defaultdict(list) for n in range(Ntimes)]
    means_gas=[defaultdict(list) for n in range(Ntimes)]
    #for nprofile, profile in enumerate(profs):
    #    print(profile)

    profiles = obj.profiles_gas
    core_list=list(profiles.keys())
    #core_list=core_list[:1] #kludge
    for core_id in core_list:
        frames = sorted(list(profiles[core_id].keys()))
        figx,axx=plt.subplots(2,4)
        for nframe,frame in enumerate(frames):

            Pgas = obj.profiles_gas[core_id][frame]
            Ppart= obj.profiles_part[core_id][frame]

            if 0:
                #density histograms
                print('kludge rho')
                continue
                rho_bins=np.geomspace(1e-3,1e8)
                PDF=True
                ax.hist(Pgas['rho_sort'], bins=rho_bins,weights=Pgas['dv_sort'],   histtype='step',density=PDF, color=[1.0,0.0,0.0],alpha=0.8)
                ax.hist(Ppart['rho_sort'], bins=rho_bins,weights=Ppart['dv_sort'], histtype='step',density=PDF, color=[0.0,1.0,0.0],alpha=0.8)
            if 1:
                PDF=True
                vr_bins=np.linspace(-20,20,32)
                #ax.hist(Pgas['vr_sort'],  bins=vr_bins,weights=Pgas['dv_sort'],  histtype='step',density=PDF, color=[1.0,0.0,0.0],alpha=0.8)
                #ax.hist(Ppart['vr_sort'], bins=vr_bins,weights=Ppart['dv_sort'], histtype='step',density=PDF, color=[0.0,1.0,0.0],alpha=0.8)
                means_par[nframe]['vr'].append(Ppart['vr_sort'].mean())
                means_gas[nframe]['vr'].append(Pgas['vr_sort'].mean())
                means_par[nframe]['vt'].append(Ppart['vt_sort'].mean())
                means_gas[nframe]['vt'].append(Pgas['vt_sort'].mean())
                #means_par[nframe]['v3d'].append(Ppart['vt_sort'].mean())
                #means_gas[nframe]['v3d'].append(Pgas['vt_sort'].mean())
                axx[0][nframe].hist(nar(Pgas['vt_sort']).flatten(),histtype='step')
                MeanVT=means_gas[nframe]['vt'][-1]
                axx[0][nframe].axvline(MeanVT)
                axx[0][nframe].set(title='%0.2f'%(MeanVT))

                axx[1][nframe].hist(nar(Ppart['vt_sort']).flatten(),histtype='step')
                MeanVT=means_par[nframe]['vt'][-1]
                axx[1][nframe].axvline(MeanVT)
                axx[1][nframe].set(title='%0.2f'%(MeanVT))
        print(means_gas[-1]['vr'][-1], MeanVT)
        outname='plots_to_sort/vt_dist_%s_c%04d.png'%(obj.this_looper.sim_name,core_id)
        figx.savefig(outname)
        plt.close(figx)
        print(outname)

    if 1:
        row = 0
        vr_bins = np.arange(-8,3,0.5)
        vt_bins = np.arange(0,8,0.5)
        for nframe in range(Ntimes):
            ax=axes[row][nframe]
            cuml=False
            ax.hist(means_par[nframe]['vt'], histtype='step',color='purple', bins=vt_bins, cumulative=cuml, density=True)
            ax.hist(nar(means_gas[nframe]['vt']).flatten(), histtype='step',color='g', bins=vt_bins, cumulative=cuml, density=True)
            print("VT",nar(means_par[nframe]['vt']).flatten())
    if 1:
        row = 1
        bins = np.arange(-8,2,0.5)
        for nframe in range(Ntimes):
            ax=axes[row][nframe]
            #qp = nar(means_par[nframe]['vr'])+0
            #qp.sort()
            #qg = nar(means_gas[nframe]['vr'])+0
            #qg.sort()
            #yp = np.arange( qp.size)/qp.size
            #yg = np.arange( qg.size)/qg.size
            #axes[1][nframe].plot( qp, yp, c='purple')
            #axes[1][nframe].plot( qg, yg, c='g')

            cuml=False
            ax.hist(means_par[nframe]['vr'], histtype='step',color='purple', bins=bins, cumulative=cuml, density=True)
            ax.hist(nar(means_gas[nframe]['vr']).flatten(), histtype='step',color='g', bins=bins, cumulative=cuml, density=True)
            #print(nar(means_gas[nframe]['vr']).flatten())

    if 1:
        row=0
        for nax, ax in enumerate(axes[row]):
            ax.set(xscale='linear',yscale='linear',xlabel='vt', ylabel='P(vt)')

    if 'rho_hist' in row_dict  and False:
        row=row_dict['rho_hist']
        for nax, ax in enumerate(axes[row]):
            ax.set(xscale='log',yscale='log',xlabel='rho', ylabel='P(rho)')
    if 'vr_hist' in row_dict:
        row=row_dict['vr_hist']
        for nax, ax in enumerate(axes[row]):
            ax.set(xscale='linear',yscale='linear',xlabel='vr', ylabel='P(vr)')
    for row in axes:
        for ax in row[1:]:
            ax.set(ylabel='', yticks=[])

    fig.savefig('plots_to_sort/sub_or_super_%s'%obj.this_looper.sim_name)
    plt.close(fig)


sim_list=['u501','u502','u503']
if 'mp' not in dir():
    mp={}
for sim in sim_list:
    if sim not in mp:
        all_cores=np.unique( TL.loops[sim].tr.core_ids)
        core_list=list(all_cores)
        core_list=None#[323]
        core_list=[323]
        core_list=[25]
        core_list=[74]
        core_list = TL.loops[sim].core_by_mode['Alone']
        #core_list=core_list[:2]
        this=radial.multipro(TL.loops[sim])
        timescale = 2 #0= 0-tsing, 1=tsing-tsing 2=4 panel
        this.run(core_list=core_list,tsing=tsing_tool[sim], timescale=timescale, get_particles=True, save_sorts=True)#, r_inflection=anne.inflection[sim])
        mp[sim]=this
if 0:
    for sim in sim_list:
        #phasor(mp[sim])
        phasor_vr_vt(mp[sim])
if 1:
    for sim in sim_list:
        souper(mp[sim])
