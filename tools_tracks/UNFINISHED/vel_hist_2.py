
from starter2 import *
import davetools as DT
reload(DT)
plt.close('all')
file_list=glob.glob('/home/dcollins/scratch/Paper19/particle_error/particle_error_test_c0031_threeframes.h5')
file_list=glob.glob('/home/dcollins/scratch/Paper19/particle_error/track_indfix_sixteenframe_core_0031.h5')
file_list=glob.glob('/home/dcollins/scratch/Paper19/particle_error/track_indfix_sixteenframe_core_0031.h5')
file_list=glob.glob('/scratch1/dcollins/Paper19/track_index_fix/track_indfix_sixteenframe_core_*.h5')
file_list=glob.glob('/home/dcollins/scratch/Paper19/track_index_fix/track_indfix_sixteenframe_core_*.h5')
if 'ext_v' not in dir():
    ext_v=DT.extents()
    ext_r=DT.extents()
    print("Running Extents")
    for nfile,fname in enumerate(file_list):
        this_looper=looper.core_looper(directory=directory)
        trw.load_loop(this_looper,fname)
        if True:
            for frame in this_looper.snaps:
                for core_id in this_looper.snaps[frame]:
                    snap = this_looper.snaps[frame][core_id]
                    if snap.R_mag.size > 1:
                        ext_v( snap.V_radial)
                        ext_r( snap.R_mag)

nvel_bins = 10
nrad_bins = 20
#velbins = np.logspace(np.log10(ext_v[0]),np.log10(ext_v[1]),nvel_bins+1)
vel_linthresh=0.1
vel_max = 50 #signed
Nlin=2 #bin edges
nlog = (nvel_bins-Nlin)//2+1
vba= np.linspace(0, vel_linthresh, Nlin)[:-1]
vbb= np.logspace(np.log10(vel_linthresh),np.log10(vel_max),nlog)
velbins0=np.concatenate([vba,vbb])
velbins = np.concatenate([-velbins0[::-1],velbins0[1:] ])
rmin = 0.05*1./128*0.5**4
rmax = 1
radbins = np.concatenate([np.zeros(1),np.logspace(np.log10(rmin),np.log10(rmax),nrad_bins)])
velhist_global = np.zeros([nvel_bins,nrad_bins])
for nfile,fname in enumerate(file_list):
    #0164.h5
    t1 = fname.split("/")[-1]
    #l = len("track_three_to_test_core_")
    #l = len("track_sixteen_frames_core_")
    l = len("track_indfix_sixteenframe_core_")

    this_cor = int(t1[l:l+4]) #[fname.index('_'):]
    #this_cor=31
    #if this_cor not in  [12]:#, 31]:
    #    continue
    print(this_cor)
    this_looper=looper.core_looper(directory=directory)
    trw.load_loop(this_looper,fname)
    thtr = this_looper.tr
    all_cores = np.unique(thtr.core_ids)
    core_list=all_cores
    rm = rainbow_map(len(all_cores))
    if 1:
        #big histogram
        asort =  np.argsort(thtr.times)
        n0=asort[0]
        tsorted = thtr.times[asort]
        velhist = np.zeros([nvel_bins,nrad_bins])

        fig=plt.figure(figsize=(4,4))
        axa=fig.subplots(1,1)
        frame_list=sorted(list(this_looper.snaps.keys()))
        tmap=rainbow_map(len(this_looper.snaps.keys()))

        for iframe,frame in enumerate(frame_list):
            for core_id in this_looper.snaps[frame]:
                snap = this_looper.snaps[frame][core_id]
                if snap.V_radial is None:
                    snap.compute_relative_coords()
                if len(snap.R_mag) < 3:
                    continue
                h, xe, ye = np.histogram2d(snap.R_mag,snap.V_radial, bins=(radbins,velbins))
                h=h.T

                rrr, vvv = np.meshgrid(xe,ye)
                velhist += h
                velhist_global += h
        if len(snap.R_mag) < 3:
            continue
        norm = mpl.colors.LogNorm(vmin=1,vmax=velhist.max())
        cmap=mpl.cm.jet
        cmap.set_under('w')
        p=axa.pcolormesh(rrr,vvv,velhist,norm=norm)
        cbar=fig.colorbar(p)

        for iframe,frame in enumerate(frame_list):
            for core_id in this_looper.snaps[frame]:
                snap = this_looper.snaps[frame][core_id]
                if len(snap.R_mag) < 3:
                    continue
                #rrr, vvv = np.meshgrid(radbins,velbins)
                #axa.scatter(snap.R_mag,snap.V_radial,c=tmap(iframe,snap.R_mag.size),s=0.1,label=str(frame))
                axa.scatter(snap.R_mag,snap.V_radial,c=tmap(iframe,snap.R_mag.size),s=0.1,label=str(frame))
        DT.axbonk(axa,xscale='log',yscale='linear',xlabel='R_mag',ylabel='V_rel',
                  xlim=ext_r, ylim=ext_v)
        axa.set_yscale('symlog',linthreshy=vel_linthresh)
        axa.set_xscale('symlog',linthreshx=2*rmin)
        outname = 'image_tracks/vel_hist_c%04d.png'%core_id
        print(outname)
        #axa.legend(loc=0)
        fig.savefig(outname)

        plt.close('all')
        fig,ax=plt.subplots(1,1)
        norm = mpl.colors.LogNorm(vmin=1,vmax=velhist_global.max())
        cmap=mpl.cm.jet
        cmap.set_under('w')
        p=ax.pcolormesh(rrr,vvv,velhist_global,norm=norm)
        DT.axbonk(ax,xscale='log',yscale='linear',xlabel='R_mag',ylabel='V_rel',
                  xlim=ext_r, ylim=ext_v)
        ax.set_yscale('symlog',linthreshy=vel_linthresh)
        ax.set_xscale('symlog',linthreshx=2*rmin)
        cbar=fig.colorbar(p)
        fig.savefig('bigtest.png')
if 0:
    if 0:
        #time plots
        asort =  np.argsort(thtr.times)
        n0=asort[0]
        tsorted = thtr.times[asort]
        nbins=30
        velhist = np.zeros([nbins, len(tsorted)])
        bins = np.logspace(np.log10(ext_v[0]),np.log10(ext_v[1]),nbins+1)

        tmap=rainbow_map(len(this_looper.snaps.keys()))
        fig=plt.figure(figsize=(4,4))
        axa=fig.subplots(1,1)
        frame_list=sorted(list(this_looper.snaps.keys()))
        for iframe,frame in enumerate(frame_list):
            for core_id in this_looper.snaps[frame]:
                snap = this_looper.snaps[frame][core_id]
                if snap.V_radial is None:
                    snap.compute_relative_coords()

                axa.scatter(snap.R_mag,snap.V_radial,c=tmap(iframe,snap.R_mag.size),s=0.1,label=str(frame))
        DT.axbonk(axa,xscale='log',yscale='linear',xlabel='R_mag',ylabel='V_rel',
                  xlim=ext_r, ylim=ext_v)
        axa.set_yscale('symlog',linthreshy=0.1)
        outname = 'image_tracks/vel_scatter_c%04d.png'%core_id
        print(outname)
        #axa.legend(loc=0)
        fig.savefig(outname)
                
