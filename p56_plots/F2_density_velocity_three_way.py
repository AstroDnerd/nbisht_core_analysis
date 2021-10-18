
from starter2 import *
import matplotlib.image as mpimg

import data_locations as dl
reload(dl)
reload(trackage)
plt.close('all')

form='pdf'

import means_etc
reload(means_etc)
#import three_loopers_mountain_top as TLM
import three_loopers_tenfour as TL4
if 'm1' not in dir():
    #m1 = means_etc.means_etc(tl.looper1 )
    #m2 = means_etc.means_etc(tl.looper2 )
    #m3 = means_etc.means_etc(tl.looper3 )
    #m1 = means_etc.means_etc(TLM.loops['u301'] )
    #m2 = means_etc.means_etc(TLM.loops['u302'] )
    #m3 = means_etc.means_etc(TLM.loops['u303'] )
    m1 = means_etc.means_etc(TL4.loops['u401'] )
    m2 = means_etc.means_etc(TL4.loops['u402'] )
    m3 = means_etc.means_etc(TL4.loops['u403'] )
    #mi = {'u301':m1,'u302':m2,'u303':m3}
    mi = {'u401':m1,'u402':m2,'u403':m3}



   
if 1:
    dbins = np.logspace( np.log10( 0.1), np.log10( 10), 16)
    vbins = np.logspace( np.log10( 1), np.log10(100),16)
    dbins = np.logspace( -0.5, 1, 10)
    vbins = np.logspace( -0.1, np.log10(20), 10)
    scales='log'
    def log_or_not(arr):
        return arr
    def log_or_notv(arr):
        return np.sqrt(arr)
else:
    scales='linear'
    dbins = None #np.linspace( ( 0.1), ( 10), 16)
    vbins = None #np.linspace( ( 1), (100),16)

    def log_or_not(arr):
        return np.log10(arr)
    def log_or_notv(arr):
        return np.log10(arr)/2

dlim = [dbins.min(), dbins.max()]
vlim = [vbins.min(), vbins.max()]
import colors
if 0:
    if 1:
        c1 = colors.color['u401']
        c2 = colors.color['u402']
        c3 = colors.color['u403']
        l1 = 'u401'
        l2 = 'u402'
        l3 = 'u403'
        ax, ax_den_hist,ax_vel_hist=means_etc.three_way_bean() 
        ok = slice(None)
        ok1 = m1.npart > 1
        ok2 = m2.npart > 1
        ok3 = m3.npart > 1
        ax.scatter(log_or_not(m1.dmeans[ok1]),log_or_notv(m1.variance[ok1]),c=c1,label=l1, s=0.1)
        ax.scatter(log_or_not(m2.dmeans[ok2]),log_or_notv(m2.variance[ok2]),c=c2,label=l2, s=0.1)
        ax.scatter(log_or_not(m3.dmeans[ok3]),log_or_notv(m3.variance[ok3]),c=c3,label=l3, s=0.1)

        r1 = ax_den_hist.hist(log_or_not(m1.dmeans[ok1]), histtype='step',color=c1,label=l1,bins=dbins)
        r2 = ax_den_hist.hist(log_or_not(m2.dmeans[ok2]), histtype='step',color=c2,label=l2,bins=dbins)
        r3 = ax_den_hist.hist(log_or_not(m3.dmeans[ok3]), histtype='step',color=c3,label=l3,bins=dbins)

        v1 = ax_vel_hist.hist(log_or_notv(m1.variance[ok1]), histtype='step', orientation='horizontal',color=c1,bins=vbins)
        v2 = ax_vel_hist.hist(log_or_notv(m2.variance[ok2]), histtype='step', orientation='horizontal',color=c2,bins=vbins)
        v3 = ax_vel_hist.hist(log_or_notv(m3.variance[ok3]), histtype='step', orientation='horizontal',color=c3,bins=vbins)
        axbonk(ax,yscale=scales, xscale=scales,  ylabel=r'$\sigma_v$', xlabel=r'$ \langle\rho\rangle$',
              xlim = dlim, ylim = vlim)
        axbonk(ax_vel_hist,yscale=scales,  xscale='linear',  ylabel=None, xlabel=r'$N$')
        axbonk(ax_den_hist,yscale='linear', xscale=scales,  ylabel=r'$N$', xlabel=None)
        ax_den_hist.legend(loc=1)
        ax_den_hist.set_xticks([])
        ax_vel_hist.set_yticks([])
        ax_den_hist.set_xlim(ax.get_xlim())
        ax_vel_hist.set_ylim(ax.get_ylim())

        odir = "plots_to_sort"
        plt.savefig(odir+'/pre_rho_rms_v_rms.%s'%form)

if 0:
    peaks = {}
    for this_simname in ['u301','u302','u303']:
        mountain_top_fname = "datasets_small/%s_mountain_tops_take_9.h5"%this_simname
        TLM.loops[this_simname].read_targets(mountain_top_fname)
        peaks[this_simname] = [TLM.loops[this_simname].targets[core_id].peak_density for core_id in mi[this_simname].core_list]
if 1:
    peaks = {}
    for this_simname in ['u401','u402','u403']:
        peaks[this_simname] = [TL4.loops[this_simname].targets[core_id].peak_density for core_id in mi[this_simname].core_list]

if 0:
    fig,ax=plt.subplots(1,1)
    for this_simname in ['u301','u302','u303']:
        c=colors.color[this_simname]
        ax.scatter( mi[this_simname].dmeans, peaks[this_simname], c=c)
    axbonk(ax,xscale='log',yscale='log',ylabel=r'$\rho_{\rm{peak}}$', xlabel=r'$\langle \rho(t=0) \rangle$')
    

    fig.savefig('plots_to_sort/dmeans.png')


if 1:
    c1 = colors.color['u401']
    c2 = colors.color['u402']
    c3 = colors.color['u403']
    l1 = 'u401'
    l2 = 'u402'
    l3 = 'u403'
    ax, ax_den_hist,ax_vel_hist=means_etc.three_way_bean()
    ok = slice(None)
    ok1 = m1.npart > 1
    ok2 = m2.npart > 1
    ok3 = m3.npart > 1
    ax.scatter(log_or_not(m1.dmeans[ok1]),log_or_not(peaks[l1]),c=c1,label=l1, s=0.1)
    ax.scatter(log_or_not(m2.dmeans[ok2]),log_or_not(peaks[l2]),c=c2,label=l2, s=0.1)
    ax.scatter(log_or_not(m3.dmeans[ok3]),log_or_not(peaks[l3]),c=c3,label=l3, s=0.1)

    r1 = ax_den_hist.hist(log_or_not(m1.dmeans[ok1]), histtype='step',color=c1,label=l1,bins=dbins)
    r2 = ax_den_hist.hist(log_or_not(m2.dmeans[ok2]), histtype='step',color=c2,label=l2,bins=dbins)
    r3 = ax_den_hist.hist(log_or_not(m3.dmeans[ok3]), histtype='step',color=c3,label=l3,bins=dbins)

    pbins = np.logspace( np.log10( 1100), np.log10( 5e7), 16)
    v1 = ax_vel_hist.hist(log_or_not(peaks[l1]), histtype='step', orientation='horizontal',color=c1,bins=pbins)
    v2 = ax_vel_hist.hist(log_or_not(peaks[l2]), histtype='step', orientation='horizontal',color=c2,bins=pbins)
    v3 = ax_vel_hist.hist(log_or_not(peaks[l3]), histtype='step', orientation='horizontal',color=c3,bins=pbins)
    axbonk(ax,yscale=scales, xscale=scales,  ylabel=r'$\rho_{peak}$', xlabel=r'$ \langle\rho\rangle$')
          #xlim = dlim, ylim = vlim)
    axbonk(ax_vel_hist,yscale=scales,  xscale='linear',  ylabel=None, xlabel=r'$N$')
    axbonk(ax_den_hist,yscale='linear', xscale=scales,  ylabel=r'$N$', xlabel=None)
#   ax_den_hist.legend(loc=1)
#   ax_den_hist.set_xticks([])
#   ax_vel_hist.set_yticks([])
#   ax_den_hist.set_xlim(ax.get_xlim())
#   ax_vel_hist.set_ylim(ax.get_ylim())

    odir = "plots_to_sort"
    plt.savefig(odir+'/pre_rho_peak_rho.%s'%form)

