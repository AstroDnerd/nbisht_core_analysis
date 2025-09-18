
from starter2 import *

import mountain_top
reload(mountain_top)
reload(looper)
import tracks_read_write
reload(tracks_read_write)
#
# get mountain tops
#

#this_simname = 'u302'
core_id = 14
#mountain_top_name = "u301_mountain_overlap_test_2.h5"
do_mountain_projections=False

def verify_cores_u302(leaf, peak_density, peak_id):
    out = True
    if peak_density < 1100:
        out = False
    return out

if 'leaf_storage' not in dir():#'leaf_storage' not in dir() and False:
    #kludge={'peak_id':[0,1,364, 113,u302_was258]})
    #kludge={'peak_id':[258]}
    #kludge={'peak_id':[core_id]}
    kludge={}
    #kludge={'peak_id':[89,11,14,350, 349, 368, 369]}
    #kludge={'peak_id':[ 368, 369]}
    leaf_storage={}
    MT=mountain_top.cut_mountain_top( this_simname, target_fname = mountain_top_name, do_projections=do_mountain_projections, 
                                     verify=verify_cores_u302,kludge=kludge, leaf_storage=leaf_storage)
                                     #density_cut_fname="u302_contour_mins.h5")
if 'peak_rho1' not in dir():

    peak_rho1=mountain_top.get_peak_densities( this_simname, target_fname = mountain_top_name, do_projections=do_mountain_projections, 
                                     verify=verify_cores_u302,kludge=kludge, leaf_storage=leaf_storage)
                                     #density_cut_fname="u302_contour_mins.h5")

if 1:
    fig,ax = plt.subplots(2,2)
    fig.subplots_adjust(wspace=0, hspace=0)
    ax0 = ax[0][0]; ax1=ax[0][1]
    ax2 = ax[1][0]; ax3=ax[1][1]
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")
    ax3.yaxis.tick_right()
    ax3.yaxis.set_label_position("right")
    ax1.xaxis.set_label_position('top') 
    ax1.xaxis.tick_top()

    pid = nar([i for i in leaf_storage])
    peak_rho = nar([np.log10( peak_rho1[i]) for i in pid])
    mass=  nar([np.log10(leaf_storage[i][1]['cell_mass'].sum().v) for i in pid])
    volume=nar([np.log10(leaf_storage[i][1]['cell_volume'].sum().v) for i in pid])
    Rmax=nar([np.log10(leaf_storage[i][1]['radius'].max().v) for i in pid])
    density=mass/volume

    ax0.xaxis.set_label_position('top') 
    ax0.xaxis.tick_top()

    ok = density < 0.5
    pfit = np.polyfit(mass[ok],density[ok],1)
    cut_line = pfit[0]*mass+pfit[1]+0.05
    ax0.plot( mass, cut_line)

    ok1 = density > cut_line
    ok2 = density <= cut_line
    ax0.scatter(mass[ok1], (mass/volume)[ok1], c='r')
    ax0.scatter(mass[ok2], (mass/volume)[ok2], c='g')
    axbonk(ax0,xlabel='mass',ylabel='density')

    ax1.scatter(mass,peak_rho)
    axbonk(ax1,xlabel='mass',ylabel='rho_max')

    ax3.scatter(Rmax,mass/volume  )
    axbonk(ax3,xlabel='Rmax',ylabel='density')


    ax2.scatter( density, peak_rho)
    axbonk(ax2,xlabel='rho',ylabel='rho_max')
    fig.savefig('plots_to_sort/%s_mass_volume.png'%this_simname)
    print('asssf')

if 0:
    fptr = open("plots_to_sort/%s_core_by_volume_below.html"%this_simname,'w')
    pid = pid[ok2]

    fptr.write("<table><tr><th> pid</th><th> Volume </th><th>img</th></td>/n")
    for nparticle,particle_id in enumerate(pid):
        fptr.write("<tr><td>%d</td><td>%f</td><td><img src='u302_peak_p%04d_34_Projection_y_density.png'</img></td></tr>\n"%(particle_id, vols[nparticle],particle_id))
    fptr.write("</table")
    fptr.close()




