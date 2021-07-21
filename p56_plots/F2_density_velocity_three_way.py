
from starter2 import *
import matplotlib.image as mpimg

import data_locations as dl
reload(dl)
reload(trackage)
plt.close('all')

form='pdf'

import means_etc
reload(means_etc)
if 'm1' not in dir():
    m1 = means_etc.means_etc(tl.looper1 )
    m2 = means_etc.means_etc(tl.looper2 )
    m3 = means_etc.means_etc(tl.looper3 )

tl.looper1.c='r'
tl.looper2.c='g'
tl.looper3.c='b'


   
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

if 1:
    if 1:
        ax, ax_den_hist,ax_vel_hist=means_etc.three_way_bean()
        ok = slice(None)
        ok1 = m1.npart > 10
        ok2 = m2.npart > 10
        ok3 = m3.npart > 10
        ax.scatter(log_or_not(m1.dmeans[ok1]),log_or_notv(m1.variance[ok1]),c=tl.looper1.c,label=tl.looper1.out_prefix, s=0.1)
        ax.scatter(log_or_not(m2.dmeans[ok2]),log_or_notv(m2.variance[ok2]),c=tl.looper2.c,label=tl.looper2.out_prefix, s=0.1)
        ax.scatter(log_or_not(m3.dmeans[ok3]),log_or_notv(m3.variance[ok3]),c=tl.looper3.c,label=tl.looper3.out_prefix, s=0.1)

        r1 = ax_den_hist.hist(log_or_not(m1.dmeans[ok1]), histtype='step',color=tl.looper1.c,label=tl.looper1.out_prefix,bins=dbins)
        r2 = ax_den_hist.hist(log_or_not(m2.dmeans[ok2]), histtype='step',color=tl.looper2.c,label=tl.looper2.out_prefix,bins=dbins)
        r3 = ax_den_hist.hist(log_or_not(m3.dmeans[ok3]), histtype='step',color=tl.looper3.c,label=tl.looper3.out_prefix,bins=dbins)

        v1 = ax_vel_hist.hist(log_or_notv(m1.variance[ok1]), histtype='step', orientation='horizontal',color=tl.looper1.c,bins=vbins)
        v2 = ax_vel_hist.hist(log_or_notv(m2.variance[ok2]), histtype='step', orientation='horizontal',color=tl.looper2.c,bins=vbins)
        v3 = ax_vel_hist.hist(log_or_notv(m3.variance[ok3]), histtype='step', orientation='horizontal',color=tl.looper3.c,bins=vbins)
        axbonk(ax,yscale=scales, xscale=scales,  ylabel=r'$\sigma_v$', xlabel=r'$ \langle\rho\rangle$',
              xlim = dlim, ylim = vlim)
        axbonk(ax_vel_hist,yscale=scales,  xscale='linear',  ylabel=None, xlabel=r'$N$')
        axbonk(ax_den_hist,yscale='linear', xscale=scales,  ylabel=r'$N$', xlabel=None)
        ax_den_hist.legend(loc=1)
        ax_den_hist.set_xticks([])
        ax_vel_hist.set_yticks([])
        ax_den_hist.set_xlim(ax.get_xlim())
        ax_vel_hist.set_ylim(ax.get_ylim())

        plt.savefig(odir+'/pre_rho_rms_v_rms.%s'%form)

