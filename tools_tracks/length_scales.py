
from starter2 import *
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
reload(dl)
plt.close('all')
color={'u05':'r','u10':'g','u11':'b'}
color.update({'u201':'r','u202':'g','u203':'b'})
import three_loopers_1tff as tl

clobber=False

#get hull_tool from convex_hull.

if 'ht1' not in dir() or clobber: 
    ht1 = hull_tool(tl.looper1)
if 'ht2' not in dir() or clobber:
    ht2 = hull_tool(tl.looper2)
    ht2.plot_2d(frames=[0])
if 'ht3' not in dir() or clobber:
    ht3 = hull_tool(tl.looper3)


if 0:
    import sf2
    reload(sf2)
    fig,ax=plt.subplots(1,1)
    sfs = {}
    pfits = {}
    Lsonic={}
    for this_looper in [tl.looper1, tl.looper2, tl.looper3]:
        name = this_looper.out_prefix
        sfs[name] = sf2.make_sf( this_looper=this_looper,frame=0, field='velocity')

        rbins,SS = sfs[name].bin_sf2(); SS/=2*np.pi

        the_x = np.log10(rbins)
        the_y = np.log10(SS)
        ok = ( np.isnan(the_x)==False)*(np.isnan(the_y)==False)
        ok = ok * (rbins > 0)*(SS>0)
        ok = ok * (rbins < 0.25)
        the_x = the_x[ok]
        the_y = the_y[ok]
        pfit = np.polyfit(the_x,the_y, 1)
        Lsonic[name] = 10**(-pfit[1]/pfit[0])
        lab = r'$p = %0.2f L_s =$'%pfit[0]
        lab += davetools.expform(Lsonic[name])
        pfits[name]=pfit
        ax.plot(rbins, SS, color[name]+"-", label = lab)
        ax.plot( rbins[ok], 1*(rbins[ok]/Lsonic[name])**pfit[0], color[name]+":")
    ax.set_title(r'$\sigma_v^2 = c_s^2 ( \ell/L_s)^p$')
    ax.legend(loc=0)
    fig.savefig('plots_to_sort/velocity_ac.pdf')


if 0:
    import p56_plots.density_AC as AC
    reload(AC)
    if 'a1' not in dir():
        a1 = AC.ac_thing('u201'); a1.plot()
        a2 = AC.ac_thing('u202'); a2.plot()
        a3 = AC.ac_thing('u203'); a3.plot()
        acs={'u201':a1,'u202':a2,'u203':a3}

if 1:
    #
    # hull lengths with density ac
    #
    fig,ax=plt.subplots(1,1)
    for nrun,ht in enumerate([ht1,ht2,ht3]):
        name = ht.this_looper.out_prefix

        #
        # Velocity AC
        #

        rbins, Vac,L  = sfs[name].bin_ac()
        ax.plot(rbins, Vac/Vac[0], c+':')
        ax.plot([L,L],[0,0.6], c+':')
        print(L)


        #
        # Density AC
        #
        ac = acs[ht.this_looper.out_prefix]
        ax.plot(ac.binned[1],ac.binned[2]/ac.binned[2][0], c=c)
        ax.plot([ac.L,ac.L],[0,0.6], c)

        #
        # Hull lengths
        #
        c=color[ ht.this_looper.out_prefix]
        hull_lengths = nar(ht.hull_volumes)**(1./3)
        vals, bins = np.histogram(hull_lengths)
        bc = 0.5*(bins[1:]+bins[:-1])
        db = (bins[1:]-bins[:-1])

        lab = r"$%s \lambda_\rho =$ %s"%(ht.this_looper.out_prefix,davetools.expform(ac.L))

        lab += r'$L_s =$%s'%davetools.expform(Lsonic[name])

        lab += r'$L_v =$%s'%davetools.expform(L)

        ax.plot(bc,vals/(vals).sum(),color=c,label=lab,linestyle="--")
        axbonk(ax,xlabel=r'$\rm{Hull\ Length}$',ylabel=r'$\rm{N}$',ylim=[0,1.])

        #rect=patches.Rectangle((0,0),ac.L,ac.ACb[0],facecolor=[0.8]*3)
        #ax.add_patch(rect)
    ax.legend(loc=0)
    fig.savefig('plots_to_sort/hull_lengths.pdf')
    plt.close('fig')
        
if 0:
    fig,ax=plt.subplots(1,1)
    #ax.hist(hull_lengths,histtype='step',color='k',normed=True)
    vals, bins = np.histogram(hull_lengths)
    bc = 0.5*(bins[1:]+bins[:-1])
    db = (bins[1:]-bins[:-1])
    ax.plot(bc,vals/vals.sum())
    ax.plot(AC.binned[1],AC.binned[2]/AC.binned[2][0],c='r')
    rect=patches.Rectangle((0,0),AC.L,AC.ACb[0],facecolor=[0.8]*3)
    ax.add_patch(rect)
    fig.savefig('plots_to_sort/%s_sizes.pdf'%this_simname)

