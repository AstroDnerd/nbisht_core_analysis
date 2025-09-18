
from starter2 import *
import annotate_particles_4
import pcolormesh_helper as pch
reload(annotate_particles_4)

def phase_movie(looper, camera=None, fields=None, 
                    core_list=None,frame_list=None, clobber=True,
                only_sphere=True):
    """
    Plots an collection of cores in a smooth manner.
    FIRST loop over frames,
        Draws boxes around the particles,  
    THEN smooths the path of the edges.
    THEN make all the plots
    """

    tr = looper.tr
    if core_list is None:
        core_list = np.unique(tr.core_ids)
    if frame_list is None:
        frame_list = looper.frame_list
    tracker_index =  [np.where(looper.tr.frames == frame)[0][0] for frame in frame_list]
    times=nar(looper.tr.times[ tracker_index] )
    all_times=looper.tr.times


    #
    #get all the miniscrubbers at once.
    #We should speed this code up.
    #

    mini_scrubbers = {}
    for core_id in core_list:
        do_velocity=True
        ms = trackage.mini_scrubber(looper.tr,core_id, do_velocity=do_velocity)
        ms.compute_ge(core_id)
        ms.compute_ke(core_id)
        mini_scrubbers[core_id]=  ms


    #
    #Loop over all cores and get the bounding box.
    #

    camera.run(core_list, frame_list, mini_scrubbers)

    for nf,frame in enumerate(frame_list):
        it=tracker_index[nf]

        # Check to see if the image was made already,
        # and skips it if it has.
        if len(core_list) == 1:
            suffix = "c%04d"%core_list[0]
        else:
            suffix = 'multi'
        outname = "%s/%s_%s_n%04d_"%(looper.plot_directory,looper.out_prefix,suffix,frame)
        got_one = False
        if not clobber:
            if len(glob.glob( "%s*"%outname)) > 0:
                got_one=True
        if got_one and not clobber:
            print("File exists, skipping")
            continue
        ds = looper.load(frame)

        left = camera.all_left[frame]
        right = camera.all_right[frame]
        center=camera.all_center[frame]
        position_dict=camera.all_positions[frame]

        #
        # main plot loop
        #
        Rmax = np.sqrt( ( (right-left)**2).max(axis=0)).max()
        sph = ds.region(center,left,right)
        ge = np.abs(sph[YT_grav_energy_2])
        ke = np.abs(sph[YT_kinetic_energy])
        xxbins=np.geomspace(5e-3,1e7,128)
        yybins=np.geomspace(5e-3,1e7,128)
        #xxbins = np.geomspace(ke.min(),ke.max(),128)
        #yybins = np.geomspace(ge[ge>0].min(),ge.max(),128)
        hist, xbins,ybins=np.histogram2d(ke[ge>0].flatten(),ge[ge>0].flatten(),bins=[xxbins,yybins])
        fig,ax=plt.subplots(1,1)
        pch.helper(hist,xbins,ybins,ax=ax)
        axbonk(ax,xscale='log',yscale='log',xlabel='KE',ylabel='GE')
        ax.plot( xxbins,xxbins,c='k')
        ax.scatter(ms.ke[:,it],np.abs(ms.ge[:,it]), edgecolor='r',s=30, facecolor='None')
        outname='plots_to_sort/phase_%s_%s_c%04d_n%04d'%(fields[0][1],fields[1][1],core_id,frame)
        fig.savefig(outname)
        print(outname)


