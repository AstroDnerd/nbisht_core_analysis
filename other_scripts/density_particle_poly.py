
from starter2 import *

import data_locations as dl
reload(dl)
import davetools as DT
from scipy.interpolate import *
reload(DT)
plt.close('all')
#file_list=glob.glob('/home/dcollins/scratch/Paper19/particle_error/particle_error_test_c0031_threeframes.h5')
#file_list=glob.glob('/home/dcollins/scratch/Paper19/particle_error/track_indfix_sixteenframe_core_0031.h5')
#file_list=glob.glob('/home/dcollins/scratch/Paper19/particle_error/track_indfix_sixteenframe_core_0031.h5')
#file_list=glob.glob('/scratch1/dcollins/Paper19/track_index_fix/track_indfix_sixteenframe_core_*.h5')
#file_list=glob.glob('/home/dcollins/scratch/Paper19/track_index_fix/track_indfix_sixteenframe_core_*.h5')

file_list=glob.glob('%s/track_indfix_sixteenframe_core_*.h5'%dl.snapshot_location)
#file_list=file_list[0:4]
plt.close('all')
G = 1620./(4*np.pi)
slope_array = []
inter_array = []
r_final_array = []
y_den_array= []
p_10_array = []
p_11_array = []
#norm = ds.arr(1,'cm')

if 'ext_v' not in dir():
    ext_d=dl.extents()
    ext_r=dl.extents()
    print("Running Extents")
    for nfile,fname in enumerate(file_list):
        this_looper=looper.core_looper(directory=dl.enzo_directory)
        trw.load_loop(this_looper,fname)
        thtr = this_looper.tr
        if True:
            for frame in this_looper.snaps:
                for core_id in this_looper.snaps[frame]:
                    density = thtr.c([core_id],'density')
                    snap = this_looper.snaps[frame][core_id]
                    if snap.R_mag.any() >1:
                        ext_r(snap.R_mag)
                        ext_d(snap.field_values['density'])
for nfile,fname in enumerate(file_list):
    #0164.h5
    print('nfile is  = %d'%nfile)
    t1 = fname.split("/")[-1]
    #l = len("track_three_to_test_core_")
    #l = len("track_sixteen_frames_core_")
    l = len("track_indfix_sixteenframe_core_")

    this_cor = int(t1[l:l+4]) #[fname.index('_'):]
    #this_cor=31
    #if this_cor not in  [12]:#, 31]:
    #    continue
    print(this_cor)
    this_looper=looper.core_looper(directory=dl.enzo_directory)
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
        density_h = density.transpose()
        density_r = np.zeros([16,181])
        all_density =[]
        all_radius = []
        #for i in range(len(density_h)):
           # density_r[i] = density_h[15-i]

        #print('the length of density 2 is %d'%len(density_2))
        #print(density)
        y = np.zeros([len(snap.R_mag)])
        y_c = np.zeros([len(snap.R_mag)])
        #all_density = np.zeros([len(tsorted),len(snap.field_values['density'])])
        #all_rad = np.zeros([len(snap.R_mag)])
        #print('length is %d'%len(snap.field_values['density']))
        #print('the length of radius array is %d'

        fig=plt.figure(figsize=(4,4))
        axa=fig.subplots(1,1)
        frame_list=sorted(list(this_looper.snaps.keys()))
        tmap=rainbow_map(len(this_looper.snaps.keys()))
        tmap_2=rainbow_map(len(tsorted))
        #for iframe in range(len(frame_list)):
            #print(iframe)

        for iframe,frame in enumerate(frame_list):
            t_good = tsorted[iframe]
            #print('the length of R_mag is: %d ' % len(snap.R_mag))
            #print(snap.R_mag)
            #print('the value of the index is = %d'%iframe)
            #print('the value of the frame is = %d'%frame)
            #print('the value of the time is = %f'%tsorted[iframe])
            #print('the length of the density is %d'%len(snap.field_values['density']))
            #print(density_r[iframe]-snap.field_values['density'])
            #for i in range(len(snap.field_values['density'])):
                #all_density[iframe][i] = snap.field_values['density'][i]
            
            
            


            for core_id in this_looper.snaps[frame]:
                snap = this_looper.snaps[frame][core_id]
                if len(snap.R_mag) < 3:
                    continue
                ave_density = snap.field_values['density'].mean(axis=0)
                ave_R = snap.R_mag.mean(axis=0)
                all_rad=snap.R_mag
                #print(snap.R_mag)
                #y = 2*(snap.R_mag)**(-2)*t_good**2
                #y_c = (snap.R_mag)**(-1.5)
                #y = (np.pi*4*G)**(-1)*(2/(snap.R_mag)**2)
                #y_c = (np.pi*4*G)**(-1)*(0.975/2)**(0.5)*(t_good**(-0.5)*snap.R_mag**(-1.5))
                #for i in range(len(snap.field_values['density'])):
                    #all_density[iframe][i] = snap.field_values['density'][i]

                the_R = snap.R_mag
                the_R[the_R<1./2048]=1./2048
                all_radius.append(the_R)
                all_density.append(snap.field_values['density'])
                #print('the length of the_R is %d'%len(snap.field_values['density']))
                ave_the_R = the_R.mean(axis=0)
                #p1 = np.polyfit(ave_the_R, ave_density,2)
                #print(p1)

                axa.scatter(the_R,snap.field_values['density'],c=tmap(iframe,snap.R_mag.size),s=0.1,label=str(frame))
                #axa.plot(the_R,np.poly1d(p1,ave_the_R),'r-')
                y = (np.pi*4*G)**(-1)*(20/(snap.R_mag)**2)
                #axa.plot( snap.R_mag, y, c= tmap_2(iframe))
                #axa.plot(ave_R,ave_density, c = 'k',marker='*')
        all_density = np.array(all_density)
        all_radius = np.array(all_radius)
        all_density= all_density.flatten()
        all_radius = all_radius.flatten()
        density_log = np.log10(all_density)
        radius_log = np.log10(all_radius)
        if len(all_radius) > 0:
            p1 = np.polyfit(radius_log,density_log,1)
            p_10_array.append(p1[0])
            p_11_array.append(p1[1])
            slope_array.append(p1[0])
            inter_array.append(10**p1[1])
            r_final_array.append(min(all_radius))
            print(p1)
        y_den = 10**p1[1]*all_radius**p1[0] 
        if len(y_den)>0:
            y_den_array.append(max(y_den))
        axa.scatter(all_radius,y_den,s=2,c = 'k',marker='*')
        #DT.axbonk(axa,xscale=('symlog',linthreshy=1/2048),yscale='log',xlabel='R_mag',ylabel='density')
        axa.set_xscale('symlog',linthreshx=1./2048)
        axa.set_yscale('log')
        axa.set_xlim(1./2048,1./10)
        axa.set_ylim(1,10**6)
        axa.set_title('the slope is %f and the y intercep is %f' %(p1[0],10**p1[1]))
        #axa.set_yscale('symlog',linthreshy=vel_linthresh)
        #axa.set_xscale('symlog',linthreshx=2*rmin)
        outname = 'image_tracks/density_part_poly_c%04d.png'%core_id
        print(outname)
        #axa.legend(loc=0)
        fig.savefig(outname)
        if len(slope_array)==160:
            plt.clf()
#----------------------- making histogram -------------------------------------------------#
            plt.hist(slope_array,bins = 2)
            outname2 = 'image_tracks/density_histo_slope_more.png'
            print(outname2)
            plt.savefig(outname2)
#------------------------- making intercept vs slope graph ---------------------------------#
            plt.clf()
            plt.scatter(slope_array,inter_array, c = 'r', marker = '*')
            plt.xscale('linear')
            plt.yscale('linear')
            plt.xlabel('slope')
            plt.ylabel('y_intercept')
            plt.ylim([min(inter_array),max(inter_array)+1])
            outname3 = 'image_tracks/slope_vs_inter_linear.png'
            plt.savefig(outname3)
            plt.clf()
            plt.scatter(slope_array,y_den_array,c='k',marker='*')
            plt.xscale('linear')
            plt.yscale('log')
            plt.xlabel('slope')
            plt.ylabel('true y_intercept')
            plt.ylim([min(y_den_array),max(y_den_array)])
            outname4 = 'image_tracks/slope_vs_y_inter.png'
            plt.savefig(outname4)
            plt.clf()
            #p_new = np.polyfit(y_den_array,y_den_max,1)
            #print("p_new is:")
            #print(p_new)
            plt.scatter(slope_array,r_final_array,c='k',marker='*')
            #plt.scatter(y_den_array,y_den_array+p_new[1],c='r',marker='*')
            plt.xlabel('slope')
            plt.ylabel('final radius')
            #plt.title('the slope is %f'%(p_new[0]))
            plt.xscale('linear')
            plt.yscale('log')
            plt.ylim=([min(r_final_array),max(r_final_array)])
            outname5 = 'image_tracks/radius_vs_slope.png'
            plt.savefig(outname5)
