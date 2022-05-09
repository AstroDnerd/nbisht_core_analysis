from starter2 import *
import matplotlib.image as mpimg

from scipy.optimize import curve_fit
import data_locations as dl
from collections import defaultdict
reload(dl)
reload(trackage)
plt.close('all')


import alpha_tools
reload(alpha_tools)


if 'do_all_plots' not in dir():
    do_all_plots = False

#import three_loopers as TL
#import three_loopers_1tff as TL
#import three_loopers_mountain_top as TLM
#import three_loopers_tenfour as TL4
import three_loopers_six as TL
import sf2
frame=0
sim_list = ['u601']

method = "vtc"
#total or radial (tr)
#central, mean, central dentisy (c,m,cd)
# vrc is radial central, for comparison with S2L.
# vtc is probably what we want for Subsonic Length and Virial Length.

if 'Lsubtool' not in dir():
    Lsubtool={}
    for this_simname in sim_list:
        Lsubtool[this_simname] = alpha_tools.sub_trial( TL.loops[this_simname])
        Lsubtool[this_simname].run(do_plots=['virial_radius'], velocity_method=method)
        #Lsubtool[this_simname].v_hist(frame, core_list = [13,14,15,16,17,18,19,21])
#   run2 = sub_trial(TLM.loops['u302'])
#   run2.v_hist(frame)#, core_list = [10,32,323])
#   run3 = sub_trial(TLM.loops['u303'])
#   run3.v_hist(frame)#, core_list = [10,32,323])
import heat_map
reload(heat_map)

if 0:
    #Lsubsonic only
    for this_simname in sim_list:
        fig,ax=plt.subplots(1,1)
        bins = np.geomspace(1e-4,0.3,32)
        qmatrix=heat_map.plot_heat( tool=Lsubtool[this_simname], quan_dict = Lsubtool[this_simname].subsonic_length, ax=ax,bins=bins)
        outname = "plots_to_sort/%s_subsonic_heatmap_%s.png"%(this_simname,method)
        axbonk(ax,xlabel=r'$t/t_{ff}$', ylabel=r'$L_{sub}$',yscale='log')
        fig.savefig(outname)
        plt.close(fig)

subvirial_ext = extents()
subsonic_ext = extents()
for core in Lsubtool[this_simname].subvirial_length:
    subvirial_ext( nar(Lsubtool[this_simname].subvirial_length[core]))
    subsonic_ext( nar(Lsubtool[this_simname].subsonic_length[core]))
if 1:
    for this_simname in sim_list:
        fig,axes=plt.subplots(2,2)
        axlist=axes.flatten()
        bins = np.concatenate([[0],np.geomspace(0.5/2048,0.6,32)])
        qmatrix=heat_map.plot_heat( tool=Lsubtool[this_simname], quan_dict = Lsubtool[this_simname].subsonic_length, 
                                   ax=axlist[0],bins=bins)
        qmatrix=heat_map.plot_heat( tool=Lsubtool[this_simname], quan_dict = Lsubtool[this_simname].subvirial_length, 
                                   ax=axlist[1],bins=bins)
        qmatrix=heat_map.plot_heat( tool=Lsubtool[this_simname], quan_dict = Lsubtool[this_simname].max_length, 
                                   ax=axlist[2],bins=bins)
        qmatrix=heat_map.plot_heat( tool=Lsubtool[this_simname], quan_dict = Lsubtool[this_simname].rms_length, 
                                   ax=axlist[3],bins=bins)
        outname = "plots_to_sort/%s_subsonic_heatmap_%s.png"%(this_simname,method)
        axbonk(axlist[0],xlabel=r'$t/t_{ff}$', ylabel=r'$L_{subsonic}$',yscale='log')
        axbonk(axlist[1],xlabel=r'$t/t_{ff}$', ylabel=r'$L_{subvirial}$',yscale='log')
        axbonk(axlist[2],xlabel=r'$t/t_{ff}$', ylabel=r'$L_{max}$',yscale='log')
        axbonk(axlist[3],xlabel=r'$t/t_{ff}$', ylabel=r'$L_{rms}$',yscale='log')

        fig.savefig(outname)
        plt.close(fig)

if 1:
    for this_simname in sim_list:
        fig,axes=plt.subplots(2,2)
        axlist=axes.flatten()
        bins = np.concatenate([[0],np.geomspace(0.5/2048,0.6,32)])
        virial_max={}
        v_m_rat = extents()
        for core_id in Lsubtool[this_simname].subvirial_length:
            virial_max[core_id] = nar(Lsubtool[this_simname].subvirial_length[core_id])/nar(Lsubtool[this_simname].max_length[core_id])
            v_m_rat(virial_max[core_id])
        qmatrix=heat_map.plot_heat( tool=Lsubtool[this_simname], quan_dict = virial_max,
                                   ax=axlist[0])
#        qmatrix=heat_map.plot_heat( tool=Lsubtool[this_simname], quan_dict = Lsubtool[this_simname].subvirial_length, 
#                                   ax=axlist[1],bins=bins)
#        qmatrix=heat_map.plot_heat( tool=Lsubtool[this_simname], quan_dict = Lsubtool[this_simname].max_length, 
#                                   ax=axlist[2],bins=bins)
#        qmatrix=heat_map.plot_heat( tool=Lsubtool[this_simname], quan_dict = Lsubtool[this_simname].rms_length, 
#                                   ax=axlist[3],bins=bins)
        outname = "plots_to_sort/%s_length_ratios_%s.png"%(this_simname,method)
        axbonk(axlist[0],xlabel=r'$t/t_{ff}$', ylabel=r'$L_{subvirial}/L_{max}$')
#        axbonk(axlist[1],xlabel=r'$t/t_{ff}$', ylabel=r'$L_{subvirial}$',yscale='log')
#        axbonk(axlist[2],xlabel=r'$t/t_{ff}$', ylabel=r'$L_{max}$',yscale='log')
#        axbonk(axlist[3],xlabel=r'$t/t_{ff}$', ylabel=r'$L_{rms}$',yscale='log')
#
        fig.savefig(outname)
        plt.close(fig)

if 0:
    fig,ax=plt.subplots(1,1)
    nsub = []
    for rrr in [run1]:
        for core_id in rrr.cores_used:
            ax.plot( rrr.times, rrr.subsonic_length[core_id],  c=[0.5]*3, linewidth=0.1)
            nsub.append( (nar(rrr.subsonic_length[core_id])>0).sum())
        axbonk(ax, xlabel='t',ylabel=r'$L_{sub}$')
        ax.set_yscale('symlog',linthresh=1./2048)
        fig.savefig('plots_to_sort/%s_subsonic_length.png'%(rrr.name))
        ax.clear()
        ax.hist(nsub)
        fig.savefig('plots_to_sort/%s_n_subsonic.png'%(rrr.name))
