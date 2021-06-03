from starter2 import *
import data_locations as dl
import davetools
reload(davetools)
reload(trackage)

plt.close('all')


#file_list=glob.glob('%s/*h5'%dl.sixteen_frame)
file_list = ['../p19_newscripts/nub_hub.h5']

#for debug purposes you may want a reduced list 
#file_list=file_list[:3]    

if 'this_looper' not in dir():
    this_looper=looper.core_looper(directory=dl.enzo_directory)
    for nfile,fname in enumerate(file_list):
        this_looper.load_loop(fname)
        print( "Reading file %d of %d"%(nfile,len(file_list)))
    thtr = this_looper.tr
    thtr.sort_time()
    all_cores = this_looper.core_list

core_list=all_cores
rm = rainbow_map(len(all_cores))
velocity_p_x = []
velocity_p_y = []
velocity_p_z = []
vel_rel_x = []
vel_rel_y = []
vel_rel_z = []
v_mag_rel = []
x_rel = []
y_rel = []
z_rel = []
R_rel = []
n_x = []
n_y = []
n_z = []
Vel_rad = []

if 0:
    if 'rho_extents' not in dir():
        rho_extents=davetools.extents()
        r_extents=davetools.extents()
        for nc,core_id in enumerate(all_cores):
            ms = trackage.mini_scrubber(thtr,core_id)
            if ms.nparticles == 1:
                continue
            density = thtr.c([core_id],'density')
            rho_extents(density)
            r_extents(ms.r)

        for nc,core_id in enumerate(core_list):
            plt.clf()
            print("why")
            density_all = []
            radius_all = []
            log_density = []
            log_radius = []
            ave_radius = []
            ave_velocity_radial = []
            ave_radius_array = []
            ave_angular_vel = []
            time_array = []

            #miniscrubber computes distance, r^2, several other quantities
            ms = trackage.mini_scrubber(thtr,core_id)
            tmap=rainbow_map(ms.ntimes)
            if ms.nparticles == 1:
               continue

            asort =  np.argsort(thtr.times)
            density = thtr.c([core_id],'density')
            if (asort != sorted(asort)).any():
                print("Warning: times not sorted.")
            n0=asort[0]
            tsorted = thtr.times[asort]

            fig, axd1=plt.subplots(1,1)

            for n_count,n_time in enumerate(asort):
                plt.clf()
                time=thtr.times[n_time]
                nx = 2048
                if time == 0:
                    continue
                #if core_id == 31:
                    #pdb.set_trace()
                x = np.floor((thtr.c([core_id],'x')*nx)[:,n_count])#or whatever the number of zones is
                cell_volume = ms.cell_volume[:,n_count]
                y = np.floor((thtr.c([core_id],'y')*nx)[:,n_count])
                z = np.floor((thtr.c([core_id],'z')*nx)[:,n_count])
                velocity_p_x = thtr.c([core_id],'velocity_x')[:,n_count]
                velocity_p_y = thtr.c([core_id],'velocity_y')[:,n_count]
                velocity_p_z = thtr.c([core_id],'velocity_z')[:,n_count]
                density_2 = thtr.c([core_id],'density')[:,n_count]
                #cell_volume = thtr.c([core_id],'cell_volume')[:,n_count]
                if len(density_2)!=0:
                    V_bulk_x = (velocity_p_x*cell_volume).sum()/(cell_volume.sum())
                    V_bulk_y = (velocity_p_y*cell_volume).sum()/(cell_volume.sum())
                    V_bulk_z = (velocity_p_z*cell_volume).sum()/(cell_volume.sum())
                    x_c = (x*cell_volume).sum()/(cell_volume.sum())
                    y_c = (y*cell_volume).sum()/(cell_volume.sum())
                    z_c = (z*cell_volume).sum()/(cell_volume.sum())
                    x_rel = x-x_c
                    y_rel = y-y_c
                    z_rel = z-z_c
                    R_rel = (x_rel**2+y_rel**2+z_rel**2)**(0.5)
                    n_x = x_rel/R_rel
                    n_y = y_rel/R_rel
                    n_z = z_rel/R_rel
                    vel_rel_x = velocity_p_x- V_bulk_x
                    vel_rel_y = velocity_p_y- V_bulk_y
                    vel_rel_z = velocity_p_z- V_bulk_z
                    V_ra = (vel_rel_x*n_x+vel_rel_y*n_y+vel_rel_z*n_z)
                    v_mag_rel = (vel_rel_x**2+vel_rel_y**2+vel_rel_z**2)**(0.5)

                    print("this is the length of R_rel")
                    print(len(R_rel))
                    print("this is the length of v_mag_rel")
                    print(len(v_mag_rel))

                c=tmap(n_count,ms.nparticles)
                this_r=ms.r[:,n_time]+0
                this_v_mag = ms.rel_vmag[:,n_time]+0
                this_V_radiant = ms.vr_rel[:,n_time]+0
                density_used = ms.density[:,n_time]+0
                cell_volume = ms.cell_volume[:,n_time]+0
                #mass_used_tot = ms.mass_total[:,n_time]+0
                mass_used= ms.mass[:,n_time]+0
                #this_angular_moment_x = ms.angular_moment_x[n_time]+0
                #this_angular_moment_y = ms.angular_moment_y[n_time]+0
                #this_angular_moment_z = ms.angular_moment_z[n_time]+0
                this_angular_momentum_rel_x  = ms.angular_momentum_rel_x[:,n_time]+0
                this_angular_momentum_rel_y = ms.angular_momentum_rel_y[:,n_time]+0
                this_angular_momentum_rel_z = ms.angular_momentum_rel_z[:,n_time]+0
                this_angular_vx = ms.angular_v_x[:,n_time]+0
                this_angular_vy = ms.angular_v_y[:,n_time]+0
                this_angular_vz = ms.angular_v_z[:,n_time]+0
                this_angular_mag_square_scatter = np.sqrt(ms.angular_v_x[:,n_time]**2+ms.angular_v_y[:,n_time]**2+ms.angular_v_z[:,n_time]**2)
                angular_vx_ave = np.prod(this_angular_vx)**(1/len(this_angular_vx))
                angular_vy_ave = np.prod(this_angular_vy)**(1/len(this_angular_vy))
                angular_vz_ave = np.prod(this_angular_vz)**(1/len(this_angular_vz))
                angular_vel_total = np.sqrt(this_angular_vx**2 + this_angular_vy**2 + this_angular_vz**2)
                p_fit = np.polyfit(np.log10(this_r),np.log10(angular_vel_total),1)
                y_fit = 10**(p_fit[1])*this_r**p_fit[0]
                angular_vel_total_ave = np.sqrt(angular_vx_ave**2 + angular_vy_ave**2 + angular_vz_ave**2)
                angular_vel_total_ave_II = np.prod(angular_vel_total)**(1/len(angular_vel_total))
                r_ave = np.prod(this_r)**(1/len(this_r))

                ave_angular_vel.append(angular_vel_total_ave_II)

            
                
                r_un = nar(sorted(np.unique(this_r)))
                v_mag_un = nar(sorted(np.unique(this_v_mag)))
                V_radiunt_un = nar(sorted(np.unique(this_V_radiant)))
                V_radiant_ave = sum(this_V_radiant*mass_used)/sum(mass_used)
                #r_ave = sum(this_r*mass_used)/sum(mass_used)
                time_array.append(r_ave)

                ave_radius_array.append(r_ave)
                ave_velocity_radial.append(V_radiant_ave)
                #--------------------- fit quantities --------------------------#
                density_all.append(density_used)
                radius_all.append(this_r)
                #plt.scatter(r_ave,angular_vel_total_ave_II,c='r',marker = '.')
                plt.scatter(this_r,angular_vel_total,c=c,label=thtr.times[n_time],s=0.1)
                plt.plot(this_r,y_fit,c='k')
                #plt.plot(time_array,ave_angular_vel,c='r')
                #plt.plot(time,angular_vel_total_ave,c=c)
                #plt.scatter(r_ave,V_radiant_ave,c = 'r',marker='*',s=200)
                #plt.plot(r_ave,V_radiant_ave, c = 'r')
                plt.xlabel('r')
                plt.ylabel('angular_velocity_total')
                plt.yscale('log')
                plt.xscale('log')
                plt.title('w = %f + %f*r'%(10**p_fit[1],p_fit[0]))
                #plt.yscale('symlog',lintresh=1e-6)
            #plt.plot(time_array,ave_angular_vel,c='r')
                outname = '%s/angular_velocity_total_sum_scatter_linear_c%04d_f%04d'%(dl.output_directory,core_id,n_time)
                plt.savefig(outname)
                print("saved "+outname)

                #axd1.plot(r_un, 100*(r_un/1e-2)**-2,c='k',linewidth=0.1)
            density_all = np.array(density_all)
            radius_all =  np.array(radius_all)
            density_all = density_all.flatten()
            radius_all =  radius_all.flatten()
            log_density = np.log10(density_all)
            log_radius = np.log10(radius_all)
            if 0:
                if len(radius_all) != 0:
                    p1 = np.polyfit(log_radius,log_density,1)
                    axd1.scatter(ave_radius_array,ave_velocity_radial[1]*(ave_radius_array[1]/ave_radius_array)**(p1[0]+2),c='g',marker = '*')


            #davetools.axbonk(axd1,xscale='log',yscale='linear',xlabel='r',ylabel=r'$\rho$',ylim = [min(ave_velocity_radial),max(ave_velocity_radial)])
                             #xlim=r_extents.minmax, ylim=rho_extents.minmax)
            #axd1.set_ylim([min(V_radiant_ave),max(V_radiant_ave)])
            #---------------- more fit quantities---------------------#
            
            outname = '%s/vel_radial_to_look_at_core_c%04d'%(dl.output_directory,core_id)
            plt.savefig(outname)
            print("saved "+outname)
