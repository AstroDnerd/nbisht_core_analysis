
from starter2 import *
import davetools as DT
reload(DT)
plt.close('all')
#data_location = '/scratch1/dcollins/Paper19/Datasets/'
import data_locations as dl
reload(dl)
#file_list=glob.glob('%s/track_indfix_sixteenframe_core_*.h5'%dl.snapshot_location)
file_list=glob.glob('/home/dcollins/scratch/Paper19/all_frames/track_all_frames_*.h5')

#This keeps track of the min and max velocity.
if 'ext_v' not in dir():
    ext_v=DT.extents()  
    ext_r=DT.extents()
    ext_v(nar([-4.29e+01, 3.26e+01]))
    ext_r(nar([2.13e-08, 1.11e+00]))
    print("Running Extents")
    if 0:
        #in case you actually need to compute the extents
        for nfile,fname in enumerate(file_list):
            this_looper=looper.core_looper(directory=dl.enzo_directory)
            trw.load_loop(this_looper,fname)
            if True:
                for frame in this_looper.snaps:
                    for core_id in this_looper.snaps[frame]:
                        snap = this_looper.snaps[frame][core_id]
                        if snap.R_mag.size > 1:
                            ext_v( snap.V_radial)
                            ext_r( snap.R_mag)

#
# Set up velocity bins.  
# Force a symetric log in the velocity.  (Ask if this isn't clear)
#

nvel_bins = 10
nrad_bins = 20
vel_max = 50 #max signed velocity
vel_linthresh=0.1 #Between += this, the plot is linear 
Nlinear=2 #Number of linear bins
nlog = (nvel_bins-Nlinear)//2+1
vba= np.linspace(0, vel_linthresh, Nlinear)[:-1] #the linear portion
vbb= np.logspace(np.log10(vel_linthresh),np.log10(vel_max),nlog)#the log portion
velbins0=np.concatenate([vba,vbb]) #put the linear and log together
velbins = np.concatenate([-velbins0[::-1],velbins0[1:] ]) #make it symmetric
rmin = 0.05*1./128*0.5**4
rmax = 1
radbins = np.concatenate([np.zeros(1),np.logspace(np.log10(rmin),np.log10(rmax),nrad_bins)])
velhist_global = np.zeros([nvel_bins,nrad_bins])

for nfile,fname in enumerate(file_list):
    
    #this is a crude way to get the core_id for debug purpose
    #t1 = fname.split("/")[-1]
    #l = len("track_indfix_sixteenframe_core_")
    #this_core = int(t1[l:l+4]) #[fname.index('_'):]
    #this_core=31
    #if this_core not in  [12]:#, 31]:
    #    continue
    #print(this_core)

    #This builds the looper object
    this_looper=looper.core_looper(directory=dl.enzo_directory)
    trw.load_loop(this_looper,fname)

    #some short hand for ease of use.
    thtr = this_looper.tr
    all_cores = np.unique(thtr.core_ids)
    core_list=all_cores
    #for making things colorful.
    #which is nice for effect, but better for keeping track of
    #what's going on.
    rm = rainbow_map(len(all_cores)) 
    tmap=rainbow_map(len(this_looper.snaps.keys()))

    frame_list=sorted(list(this_looper.snaps.keys()))
    if 1:
        #big histogram

        #the histogram t fill
        velhist = np.zeros([nvel_bins,nrad_bins])

        #the plot tools
        fig=plt.figure(figsize=(4,4))
        axa=fig.subplots(1,1)

        for iframe,frame in enumerate(frame_list):
            for core_id in this_looper.snaps[frame]:
                snap = this_looper.snaps[frame][core_id]
                if snap.V_radial is None:
                    snap.compute_relative_coords()
                if len(snap.R_mag) < 3:
                    continue
                if snap.time == 0:
                    print("NO ZERO")
                    continue

                #makes a 2d histogram, like phase diagrams in yt.
                h, xe, ye = np.histogram2d(snap.R_mag,snap.V_radial, bins=(radbins,velbins))
                h=h.T #the transpose so that data structures line up.

                #this is the x and y coordinates of the 2d histogram that we'll
                #plot
                rrr, vvv = np.meshgrid(xe,ye)

                #For this code, we'll only plot the cumulative of all 
                #histograms.
                velhist += h
                velhist_global += h
        if len(snap.R_mag) < 3:
            continue

        #plot the histogram.
        #this is cumbersome, make sure each of these pieces
        #makes sense.
        norm = mpl.colors.LogNorm(vmin=1,vmax=velhist.max())
        cmap=mpl.cm.jet
        cmap.set_under('w')
        p=axa.pcolormesh(rrr,vvv,velhist,norm=norm)
        cbar=fig.colorbar(p)

        #plot a scatter plot along with the histogram
        for iframe,frame in enumerate(frame_list):
            for core_id in this_looper.snaps[frame]:
                snap = this_looper.snaps[frame][core_id]
                if len(snap.R_mag) < 3:
                    continue
                #rrr, vvv = np.meshgrid(radbins,velbins)
                axa.scatter(snap.R_mag,snap.V_radial,c=tmap(iframe,snap.R_mag.size),s=0.1,label=str(frame))

        #this is a dumb bit of code that shortens plot formatting
        DT.axbonk(axa,xscale='log',yscale='linear',xlabel='R_mag',ylabel='V_rel',
                  xlim=ext_r, ylim=ext_v)
        #though the dumb bit of code doesn't deal with symlog
        axa.set_yscale('symlog',linthreshy=vel_linthresh)
        axa.set_xscale('symlog',linthreshx=2*rmin)
        outname = 'image_tracks/vel_hist_c%04d.png'%core_id
        print(outname)
        fig.savefig(outname)
        print("saved "+outname)
        plt.close('all')
        
        #now we'll plot the multi-core version
        fig,ax=plt.subplots(1,1)
        norm = mpl.colors.Normalize(vmin=1,vmax=velhist_global.max())
        cmap=mpl.cm.jet
        cmap.set_under('w')
        p=ax.pcolormesh(rrr,vvv,velhist_global,norm=norm)
        DT.axbonk(ax,xscale='log',yscale='linear',xlabel='R_mag',ylabel='V_rel',
                  xlim=ext_r, ylim=ext_v)
        ax.set_yscale('symlog',linthreshy=vel_linthresh)
        ax.set_xscale('symlog',linthreshx=2*rmin)
        cbar=fig.colorbar(p)
        outname = 'bigtest.png'
        fig.savefig(outname)
        print("saved "+outname)
