
from starter2 import *
import matplotlib.image as mpimg

import matplotlib.patches as patches
import data_locations as dl
reload(dl)
plt.close('all')

#file_list=glob.glob('%s/*h5'%dl.sixteen_frame)
file_list=glob.glob("../Datasets/u10_primitives_cXXXX_n0000.h5")
this_simname = 'u10'
out_prefix = this_simname

#file_list = file_list[:2]

#import tools_tracks.alpha_properties as ap 
#import p56_plots.all_pdf as ap
#this_looper=ap.this_looper
#thtr=this_looper.tr
#thtr.sort_time()


#
# dear future self:
#    calling load_loop on the u10 dataset seems to choke.
#    Reading directly works fine. 
#
if 'this_looper' not in dir():
    this_looper=looper.core_looper(directory=dl.sims[this_simname],savefile=file_list[0])
    #for nfile,fname in enumerate(file_list):
    #    this_looper.load_loop(fname)
    #    print( "File %d of %d"%(nfile,len(file_list)))
    #thtr = this_looper.tr
    #thtr.sort_time()
all_cores = np.unique(thtr.core_ids)
rm = rainbow_map(len(all_cores))
core_list = all_cores
#core_list=[7]
fig_many, ax_many = plt.subplots(1,1)
if 'cell_volumes' not in dir():
    cell_volumes = []
    hull_volumes = []
    good_cores=[]
    for core_id in core_list:
        ms = trackage.mini_scrubber(thtr,core_id)
        if ms.r.shape[0] <= 4:
            continue
        good_cores.append(core_id)
        tmap=rainbow_map(ms.ntimes)
        asort =  np.argsort(thtr.times)
        delta=0.1

        rrr = ms.r[:,asort]
        dv = rrr[:,1:] - rrr[:,:-1]
        mask = slice(None)# dv[:,1]*dv[:,2] < 0
        for it,nt in enumerate([0]):#asort):
            mask = slice(None)


            ax_many.clear()
            ax_many.plot([0,1,1,0,0],[0,0,1,1,0])
            this_frame = sorted(thtr.frames)[nt]
            cen=this_looper.snaps[this_frame][core_id].R_centroid
            cenx = cen[0]; ceny=cen[1]
            #ax_many.scatter(ms.this_x[mask,nt],ms.this_y[mask,nt],s=0.1)
            #this_x,this_y,this_z=ms.raw_x[mask,nt],ms.raw_y[mask,nt], ms.raw_z[mask,nt]
            this_x,this_y,this_z=ms.this_x[mask,nt],ms.this_y[mask,nt], ms.this_z[mask,nt]
            ax_many.scatter(this_x, this_y ,s=0.1)
            from scipy.spatial import ConvexHull
            points_3d = np.array(list(zip(this_x,this_y,this_z)))
            hull_3d = ConvexHull(points_3d)
            hull_volumes.append(hull_3d.volume)
            cell_volumes.append( thtr.c([core_id],'cell_volume')[mask,nt].sum())
            ax_many.plot(points_3d[hull_3d.vertices,0], points_3d[hull_3d.vertices,1], 'r', lw=2)

            points_2d = np.array(list(zip(this_x,this_y)))
            hull_2d = ConvexHull(points_2d)
            ax_many.plot(points_2d[hull_2d.vertices,0], points_2d[hull_2d.vertices,1], 'k')
            #ax_many.plot(points[hull.vertices[0],0], points[hull.vertices[0],1], 'ro')
            #ax_many.scatter( cenx,ceny,c='k',s=3)
            #ax_many.set_xlim(-delta,1+delta); ax_many.set_ylim(-delta,1+delta)
            title='t=%0.5f'%thtr.times[nt]
            ax_many.set_title(title)
            ax_many.set_aspect('equal')
            outname = '%s/%s_hull_3d_t_c%04d_n%04d.png'%(dl.output_directory,this_simname,core_id,it)
            fig_many.savefig(outname)
            print("Wrote "+outname)

if 1:
    good_cores=nar(good_cores)
    hull_volumes=nar(hull_volumes)
    cell_volumes=nar(cell_volumes)
    fig,ax=plt.subplots(1,1)
    hull_lengths = nar(hull_volumes)**(1./3)
    ax.hist(hull_lengths,histtype='step',color='k')
    axbonk(ax,xlabel=r'$\rm{Hull\ Length}$',ylabel=r'$\rm{N}$')
    fig.savefig('plots_to_sort/%s_hull_lengths.pdf'%this_simname)
    
    ax.clear()
    odd = hull_volumes < cell_volumes
    not_odd = hull_volumes >= cell_volumes
    ax.scatter(hull_volumes,cell_volumes,c='k')
    ax.scatter(hull_volumes[odd],cell_volumes[odd],c='r')
    ax.set_aspect('equal')
    #axbonk(ax,xlabel=r'$\rm{Hull\ Volume}$',ylabel=r'$\rm{Cell\ Volume}$',xlim=[0,0.07],ylim=[0,0.07])
    ext=extents()
    ext(hull_volumes)
    ext(cell_volumes)
    axbonk(ax,xlabel=r'$\rm{Hull\ Volume}$',ylabel=r'$\rm{Cell\ Volume}$',xlim=ext.minmax,ylim=ext.minmax,
          xscale='log',yscale='log')
    fig.savefig('plots_to_sort/%s_volume_fractions.pdf'%this_simname)

if 1:
    import p56_plots.density_AC as AC
    #reload(AC)
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


