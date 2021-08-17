
from starter2 import *
import davetools
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
reload(dl)
plt.close('all')
import colors
reload(colors)

#import three_loopers_1tff as tl
import three_loopers_mountain_top as TLM

clobber=False

import convex_hull_tools as CHT
#get hull_tool from convex_hull.

#if 'ht1' not in dir() or clobber: 
#    ht1 = hull_tool(tl.looper1)
#if 'ht2' not in dir() or clobber:
#    ht2 = hull_tool(tl.looper2)
#    ht2.plot_2d(frames=[0])
#if 'ht3' not in dir() or clobber:
#    ht3 = hull_tool(tl.looper3)
if 'ht1' not in dir() or clobber: 
    ht1 = CHT.hull_tool(TLM.loops['u301'])
if 'ht2' not in dir() or clobber:
    ht2 = CHT.hull_tool(TLM.loops['u302'])
    #ht2.plot_2d(frames=[0])
if 'ht3' not in dir() or clobber:
    ht3 = CHT.hull_tool(TLM.loops['u303'])

looper_list = [TLM.loops['u301'], TLM.loops['u302'], TLM.loops['u303']]
ht_list = [ht1,ht2,ht3]

for ht in ht_list:
    if len(ht.hull_volumes) == 0:
        ht.make_hulls()


if 'Lsonic' not in dir():
    import sf2
    reload(sf2)
    fig,ax=plt.subplots(1,1)
    sfs = {}
    pfits = {}
    Lsonic={}
    for this_looper in looper_list:
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
        ax.plot(rbins, SS, colors.color[name]+"-", label = lab)
        ax.plot( rbins[ok], 1*(rbins[ok]/Lsonic[name])**pfit[0], colors.color[name]+":")
    ax.set_title(r'$\sigma_v^2 = c_s^2 ( \ell/L_s)^p$')
    ax.legend(loc=0)
    fig.savefig('plots_to_sort/velocity_ac.pdf')


if 1:
    import p56_plots.density_AC as AC
    reload(AC)
    if 'a1' not in dir():
        a1 = AC.ac_thing('u301'); a1.plot()
        a2 = AC.ac_thing('u302'); a2.plot()
        a3 = AC.ac_thing('u303'); a3.plot()
        acs={'u301':a1,'u302':a2,'u303':a3}

if 1:
    #
    # hull lengths with density ac
    #
    L_vel={}
    L_rho={}
    L_jea={}
    fig,ax=plt.subplots(1,1)
    for nrun,ht in enumerate([ht1,ht2,ht3]):
        c=colors.color[ ht.this_looper.out_prefix]
        name = ht.this_looper.out_prefix

        do_label=False
        lab=None
        if nrun==0:
            do_label=True

        #
        # Hull lengths
        #
        c=colors.color[ ht.this_looper.out_prefix]
        hull_lengths = nar(ht.hull_volumes)**(1./3)
        vals, bins = np.histogram(hull_lengths)
        bc = 0.5*(bins[1:]+bins[:-1])
        db = (bins[1:]-bins[:-1])

        stuff = r"$%s \lambda_\rho =$ %s"%(ht.this_looper.out_prefix,davetools.expform(ac.L))

        stuff += r'$L_s =$%s'%davetools.expform(Lsonic[name])

        stuff += r'$L_v =$%s'%davetools.expform(L)

        stuff += r'$L_J =$%s'%davetools.expform(Ljeans)

        if do_label:
            lab='Hull Lengths'
        ax.plot(bc,vals/(vals).sum(),color=c,label=lab)

        #
        # Velocity AC
        #

        if do_label:
            lab = 'velocity AC'
        rbins, Vac,L  = sfs[name].bin_ac()
        ax.plot(rbins, Vac/Vac[0], c=c,linestyle='--',label=lab)
        L_vel[name] = L
        print(L)


        #
        # Density AC
        #
        if do_label:
            lab = 'density AC'
        ac = acs[ht.this_looper.out_prefix]
        ax.plot(ac.binned[1],ac.binned[2]/ac.binned[2][0], c=c,linestyle=':',label=lab)
        #ax.plot([ac.L,ac.L],[0,0.6], c)
        L_rho[name]=ac.L

        #
        # Jeans Length
        #
        csound = 1
        rho0 = 1
        G = 1620/(4*np.pi)
        Ljeans = csound/np.sqrt(G*rho0)
        #ax.plot([Ljeans,Ljeans],[0,0.6], 'k')
        L_jea[name] = Ljeans


        #rect=patches.Rectangle((0,0),ac.L,ac.ACb[0],facecolor=[0.8]*3)
        #ax.add_patch(rect)

    #plot velocity autocorrelation length
    mean_Lvel = np.mean(list(L_vel.values()))
    ax.plot( [mean_Lvel]*2,[0,0.35],c=[0.5]*3,linestyle='--',label=r'$L_{\rm{AC,v}}$')

    #plot density autocorrelation length
    mean_Lrho = np.mean(list(L_rho.values()))
    ax.plot( [mean_Lrho]*2,[0,0.35],c=[0.5]*3,linestyle=':', label=r'$L_{\rm{AC,\rho}}$')
    #plot Jeans
    mean_Jea = np.mean(list(L_jea.values()))
    ax.scatter([mean_Jea],0.4,marker='*',c=[[0.5]*3], label=r'$L_{\rm{Jeans}}$')
    #plot Sonic Length
    mean_Sonic = np.mean(list(Lsonic.values()))
    ax.scatter([mean_Sonic],0.4,marker='*',c='k', label=r'$L_{\rm{sonic}}$')


    axbonk(ax,xlabel=r'$\rm{Hull\ Length}$',ylabel=r'$\rm{N}$',ylim=[0,1.], xlim=[0,0.55])
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

