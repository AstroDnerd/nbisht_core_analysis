"""
This just pulls partilce data from enzo and stores it.
Changing
core_list 
frame_list
fields
changes what gets extracted.
"""
from starter2 import *
import xtra_energy
import data_locations as dl
reload(dl)
#
# set sim
#
if 'this_simname' not in dir():
    this_simname = 'u11'

def three_way_bean():
    from matplotlib.ticker import NullFormatter
    nullfmt = NullFormatter()         # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(1, figsize=(8, 8))

    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    return axScatter,axHistx, axHisty

import three_loopers_1tff as tl

if 0:
    #
    # make first_last object
    #
    core_list =  [0]
    frame_list = [0, dl.target_frames[this_simname]]

    #frame_list = list(range(0, dl.target_frames[this_simname],10))+[dl.target_frames[this_simname]]
    fields = ['x','y','z','density','velocity_magnitude','magnetic_field_strength', 'velocity_divergence']
    fields += ['velocity_x','velocity_y','velocity_z']
    fields += ['magnetic_field_%s'%s for s in 'xyz']
    fields += ['PotentialField']
    output_base = "primitive_test"
    derived=[]

    output_name = '%s_first_last_t2_nXXX0.h5'%(this_simname)
    this_looper = looper.core_looper(directory= dl.sims[this_simname],
                                     sim_name = this_simname,
                                     out_prefix = this_simname,
                                     target_frame = dl.target_frames[this_simname],
                                     frame_list = frame_list,
                                     core_list =   core_list,
                                     fields_from_grid=fields,
                                     derived = derived
                                  )
    ds = this_looper.load(frame_list[0])
    ad = ds.all_data() #ds.region([0.5]*3,[0.4]*3,[0.6]*3)
    this_looper.target_indices[0]=ad['particle_index']
    this_looper.get_tracks()
    this_looper.save(output_name)

if 1:
    savefile = '%s_first_last_t2_nXXX0.h5'%(this_simname)
    new_loop = looper.core_looper(directory= dl.sims[this_simname])
    new_loop.load_loop(savefile)

    field = 'magnetic_field_strength'
    plotdir = "./plots_to_sort/"
    outname = "%s/%s_%s_first_last.pdf"%(plotdir,field,this_simname)

    rho0 = new_loop.tr.track_dict[field][:,0]
    rho1 = new_loop.tr.track_dict[field][:,-1]

    nx=ny=64
    xbins = np.logspace(np.log10(rho0.min()), np.log10(rho0.max()),nx+1)
    ybins = np.logspace(np.log10(rho1.min()), np.log10(rho1.max()),nx+1)

    this_hist, xedge, yedge= np.histogram2d(rho0, rho1, bins=[xbins,ybins])

    def cen(arr):
        return 0.5*(arr[1:]+arr[:-1])
    TheX = np.r_[(ny)*[cen(xbins)]].transpose()
    TheY = np.r_[(nx)*[cen(ybins)]]
    x_del = (xbins[1:]-xbins[:-1])
    y_del = (ybins[1:]-ybins[:-1])
    x_del.shape = (x_del.size,1)
    dv = 1./(x_del*y_del)
    pdf = this_hist/(this_hist*dv).sum()

    #fig,ax=plt.subplots(1,1)
    plt.clf()
    ax, ax_x_hist,ax_y_hist=three_way_bean()
    cmap = copy.copy(mpl.cm.get_cmap("viridis"))
    cmap.set_under('w')
    minmin = pdf[pdf>0].min()
    norm = mpl.colors.LogNorm(vmin=minmin,vmax=pdf.max())
    ploot=ax.pcolormesh(TheX, TheY, pdf,cmap=cmap,norm=norm)#,alpha=0.2)
    ax.plot([1e-3,1e2],[1e-3,1e2],c='k')

    track_array= tl.looper2.tr.track_dict[field]
    xxx = track_array[:,0]
    yyy = track_array[:,-1]
    ax.scatter(xxx,yyy,s=0.1,c='r')



    ax_x_hist.hist(rho0, histtype='step',bins=xbins)
    ax_y_hist.hist(rho1, histtype='step',bins=ybins, orientation='horizontal')
    axbonk(ax,xscale='log',yscale='log', xlabel=r'$\rho_0$', ylabel= r'$\rho_{final}$', xlim=[xbins.min(),xbins.max()],ylim=[ybins.min(),ybins.max()])
    axbonk(ax_x_hist,xscale='log',yscale='log', xlim=[xbins.min(),xbins.max()], xlabel='', ylabel=r'$N$')
    axbonk(ax_y_hist,xscale='log',yscale='log', ylim=[ybins.min(),ybins.max()], ylabel=r'$N$')
    ax_x_hist.set_xticks([])
    ax_y_hist.set_yticks([])
    ax_y_hist.yaxis.tick_right()
    plt.savefig(outname)
    print(outname)

