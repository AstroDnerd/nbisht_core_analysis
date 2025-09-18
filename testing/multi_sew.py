
from starter2 import *
import matplotlib.image as mpimg

import davetools as DT
import data_locations as dl
reload(dl)

#file_list=glob.glob('%s/track_indfix_sixteenframe_core_*.h5'%dl.snapshot_location)
gloob='/home/dcollins/scratch/Paper19/all_frames/track_all_frames_c0079.h5'
file_list=glob.glob(gloob)
plt.close('all')

#This keeps track of the min and max velocity.
if 'ext_v' not in dir():
    ext_v=DT.extents()  
    ext_r=DT.extents()
    ext_rho=DT.extents()
    ext_v(nar([-4.29e+01, 3.26e+01]))
    ext_r(nar([2.13e-08, 1.11e+00]))
    ext_rho(nar([7.95e-03, 5.43e+07]))
    print("Running Extents")
    if 0:
        #in case you actually need to compute the extents
        for nfile,fname in enumerate(file_list):
            print(fname)
            this_looper=looper.core_looper(directory=dl.enzo_directory)
            trw.load_loop(this_looper,fname)
            ext_rho(this_looper.tr.track_dict['density'])
            continue
            if True:
                for frame in this_looper.snaps:
                    for core_id in this_looper.snaps[frame]:
                        snap = this_looper.snaps[frame][core_id]
                        if snap.R_mag.size > 1:
                            ext_v( snap.V_radial)
                            ext_r( snap.R_mag)


for nfile,fname in enumerate(file_list):

    #this is crude
    #t1 = fname.split("/")[-1]
    #l = len("track_indfix_sixteenframe_core_")
    #this_cor = int(t1[l:l+4]) #[fname.index('_'):]
    #if this_cor not in  [12]:#, 31]:
    #    continue
    #print(this_cor)

    plt.clf()
    this_looper=looper.core_looper(directory=dl.enzo_directory)
    trw.load_loop(this_looper,fname)
    thtr = this_looper.tr
    all_cores = np.unique(thtr.core_ids)
    core_list=all_cores
    rm = rainbow_map(len(all_cores))

    if 1:
        #time plots
        asort =  np.argsort(thtr.times)
        n0=asort[0]
        tsorted = thtr.times[asort]
        xxx = 1./128/np.logspace(0,4,5,base=2)
        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id)
            density = thtr.c([core_id],'density')
            tmap=rainbow_map(ms.ntimes)
            for it,nt in enumerate(asort):
                frame = thtr.frames[nt]
                basedir = "/home/dcollins/RESEARCH2/Paper19_47_overlap/0000_density_tracks/x/ALL"
                basedir = "/home/dcollins/RESEARCH2/Paper19_47_overlap/0000_density_tracks/AllFrames/"
                core_dir = "%s/c%04d"%(basedir,core_id)
                #image = "%s/test_full_particles_c%04d_n%04d_Projection_x_density.png"%(core_dir,core_id,frame)
                image = "%s/test_c%04d_n%04d_Projection_x_density.png"%(core_dir,core_id,frame)
                
                image = "%s/test_c%04d_n%04d_centered_Projection_x_density.png"%(core_dir,core_id,frame)
                name_base = "image_tracks/rho_r/"
                rhoname = name_base+'/rho_r_step_c%04d_n%04d.png'%(core_id,frame)
                name_base = "image_tracks/rho_time/"
                outname = name_base+'/rho_t_fit2_c%04d_n%04d.png'%(core_id,it)#frame)
                skip = False

                for fff in [image, rhoname, outname]:
                    if not os.path.exists(fff):
                        skip=True
                        print("MISSING",fff)
                if skip:
                    #print("MISSING")
                    continue
                print("to read %s"%image)

                img1 = mpimg.imread(image) 
                img2 = mpimg.imread(outname) 
                img3 = mpimg.imread(rhoname)
                nx =  max([img1.shape[0],img2.shape[0],img3.shape[0]])
                nx = max([nx, img2.shape[0]+img3.shape[0]])
                ny = img1.shape[1]+max([img2.shape[1] , img3.shape[1]])
                #ny = img1.shape[1]+max([img2.shape[1] ])#, img2.shape[1]])
                both = np.zeros([nx,ny,4])
                both[ 0:img1.shape[0] , 0:img1.shape[1]]=img1
                both[ 0:img2.shape[0] , img1.shape[1]:(img1.shape[1]+img2.shape[1])]=img2
                s30 = img3.shape[0]
                s31 = img3.shape[1]
                x1=img2.shape[0]; x2=(img2.shape[0]+s30)#img3.shape[0])
                y1=img1.shape[1]; y2=(img1.shape[1]+s31)#img3.shape[1])
                both[ x1:x2, y1:y2 ]= img3[0:s30,0:s31]
                oname2 = "image_tracks/density_3_c%04d_n%04d"%(core_id,frame)
                print(oname2)
                plt.imsave(oname2,both)

            if 0:#for it,nt in enumerate(asort):
                oname2 = "image_tracks/density_3_c%04d_n%04d"%(core_id,frame)
                if os.path.exists(oname2):
                    continue
                plt.scatter(ms.r[:,nt], density[:,nt],c=[tmap(it)]*ms.r.shape[0],s=0.1,
                           label = "%0.4f"%thtr.times[nt])
                if it==0:
                    powerline(plt,1e-3,1e-1,5e4,-2,c='k')
                plt.xscale('symlog',linthreshx=xxx.min()/0.5);
                plt.title("%0.4f"%thtr.times[nt])
                #plt.xlim(1e-5,0.5)
                plt.xlim(0,0.5)
                plt.yscale('log')
                plt.ylim(0.01,4e6)
                #plt.legend(loc=1)
                frame = thtr.frames[nt]
                outname = 'image_tracks/rho_r_step_c%04d_n%04d.png'%(core_id,frame)
                plt.savefig(outname)

                print('saved '+outname)
