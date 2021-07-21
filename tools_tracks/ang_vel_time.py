from starter2 import *
import matplotlib.image as mpimg

from scipy.optimize import curve_fit
import data_locations as dl
from collections import defaultdict
reload(dl)
plt.close('all')



class ang_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=nar([])
        self.slopes=defaultdict(list)
        self.times=[]
    def get_extents(self,core_list=None):



        self.ext_Omega = extents()
        self.ext_r = extents()
        self.ext_I = extents()
        self.ext_j = extents()


        self.ext_Omega(nar([1,5e4]))
        self.ext_j(nar([1e-11,0.25e-4]))
        self.ext_I(nar([1e-14,2e-6]))
        return


        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores
        rm = rainbow_map(len(all_cores))

        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)

            if ms.nparticles < 3:
                continue


            j = np.sqrt( ms.angular_momentum_rel_x**2 + ms.angular_momentum_rel_y**2 + ms.angular_momentum_rel_z**2)
            self.ext_j(j)

            Omega = np.sqrt(ms.angular_v_x**2 + ms.angular_v_y**2 + ms.angular_v_z**2)
            self.ext_Omega(Omega)

            self.ext_I( ms.I_ii)






    def run(self,do_all_plots=True,core_list=None,plot_each_frame=False):

        self.get_extents(core_list)

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

        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores
        rm = rainbow_map(len(all_cores))

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        self.times=thtr.times[ thtr.times > 0] 
        for core_id in core_list:
            plt.clf()
            density_all = []
            radius_all = []
            log_density = []
            log_radius = []
            ave_radius = []
            ave_velocity_radial = []
            ave_radius_array = []
            ave_angular_vel = []
            time_array = []
            ms = trackage.mini_scrubber(thtr,core_id)

            if ms.nparticles < 3:
                print("Particle count too small ",ms.nparticles)
                continue
            print('go ', core_id)
            self.cores_used=np.append(self.cores_used,core_id)
            n0=0
            ms = trackage.mini_scrubber(thtr,core_id)
            tmap=rainbow_map(ms.ntimes)
            if ms.nparticles == 1:
               continue

            prefix = self.this_looper.out_prefix
            asort =  np.argsort(thtr.times)
            density = thtr.c([core_id],'density')
            if (asort != sorted(asort)).any():
                print("Warning: times not sorted.")
            n0=asort[0]
            tsorted = thtr.times[asort]

            fig, axd1=plt.subplots(2,2)

            for n_count,n_time in enumerate(asort):
                if plot_each_frame:
                    plt.clf()
                time=thtr.times[n_time]
                frame=thtr.frames[n_time]
                nx = 2048
                if time == 0:
                    continue
                #if core_id == 31:
                    #pdb.set_trace()
                x = np.floor((thtr.c([core_id],'x')*nx)[:,n_count])#or whatever the number of zones is
                y = np.floor((thtr.c([core_id],'y')*nx)[:,n_count])
                z = np.floor((thtr.c([core_id],'z')*nx)[:,n_count])
                cell_volume = ms.cell_volume[:,n_count]
                density_2 = thtr.c([core_id],'density')[:,n_count]
                #cell_volume = thtr.c([core_id],'cell_volume')[:,n_count]
                c=tmap(n_count,ms.nparticles)
                this_r=ms.r[:,n_time]+0
                this_r[ this_r < 0.25/2048] = 0.25/2048
                this_v_mag = ms.rel_vmag[:,n_time]+0
                this_V_radiant = ms.vr_rel[:,n_time]+0
                this_density = ms.density[:,n_time]+0
                cell_volume = ms.cell_volume[:,n_time]+0
                #mass_used_tot = ms.mass_total[:,n_time]+0
                #this_angular_moment_x = ms.angular_moment_x[n_time]+0
                #this_angular_moment_y = ms.angular_moment_y[n_time]+0
                #this_angular_moment_z = ms.angular_moment_z[n_time]+0
                this_angular_momentum_rel_x  = ms.angular_momentum_rel_x[:,n_time]+0
                this_angular_momentum_rel_y = ms.angular_momentum_rel_y[:,n_time]+0
                this_angular_momentum_rel_z = ms.angular_momentum_rel_z[:,n_time]+0
                this_j = np.sqrt( this_angular_momentum_rel_x**2+ this_angular_momentum_rel_y**2+ this_angular_momentum_rel_z**2)

                #matches this_j almost exactly.
                other_j = np.sqrt( ms.angular_moment_x[:,n_time]**2 + ms.angular_moment_y[:,n_time]**2 + \
                                  ms.angular_moment_z[:,n_time]**2 )


                this_angular_vx = ms.angular_v_x[:,n_time]+0
                this_angular_vy = ms.angular_v_y[:,n_time]+0
                this_angular_vz = ms.angular_v_z[:,n_time]+0
                #this_angular_mag_square_scatter = np.sqrt(ms.angular_v_x[:,n_time]**2+ms.angular_v_y[:,n_time]**2+ms.angular_v_z[:,n_time]**2)
                #angular_vx_ave = np.prod(this_angular_vx)**(1/len(this_angular_vx))
                #angular_vy_ave = np.prod(this_angular_vy)**(1/len(this_angular_vy))
                #angular_vz_ave = np.prod(this_angular_vz)**(1/len(this_angular_vz))
                #angular_vel_total_ave = np.sqrt(angular_vx_ave**2 + angular_vy_ave**2 + angular_vz_ave**2)
                angular_vel_total = np.sqrt(this_angular_vx**2 + this_angular_vy**2 + this_angular_vz**2)


                other_j = ms.I_ii[:,n_time]*angular_vel_total
                if this_r.sum() < 1e-16:
                    self.slopes[core_id].append(2)
                    p_fit=[2,2]
                else:
                    p_fit = np.polyfit(np.log10(this_r),np.log10(angular_vel_total),1)
                    y_fit = 10**(p_fit[1])*this_r**p_fit[0]
                    angular_vel_total_ave_II = np.prod(angular_vel_total)**(1/len(angular_vel_total))

                    ave_angular_vel.append(angular_vel_total_ave_II)

                    self.slopes[core_id].append(p_fit[0])

            
                if do_all_plots:                

                    #--------------------- fit quantities --------------------------#

                    axd1[0][0].scatter(this_r,angular_vel_total,c=c,label=thtr.times[n_time],s=0.1)
                    #axd1[0][0].plot(this_r,y_fit,c='k')
                    axd1[1][0].scatter(this_r, ms.I_ii[:,n_time], c=c, s=0.1)
                    axd1[0][1].scatter(this_r, this_j, c=c, s=0.1)
                    third_j = ms.mass[:,n_time]*this_r**2*angular_vel_total
                    
                    axd1[1][1].scatter(this_r,angular_vel_total**2/this_density,c=c,s=0.1)
                    if plot_each_frame:
                        plt.xlabel('r')
                        plt.ylabel('angular_velocity_total')
                        plt.yscale('log')
                        plt.xscale('log')
                        plt.title('w = %f + %f*r'%(10**p_fit[1],p_fit[0]))
                        outname = '%s/%s_angular_velocity_time_c%04d_n%04d'%(dl.output_directory,prefix,core_id,frame)
                        plt.savefig(outname)
                        print("saved "+outname)

            
            axd1[1][0].plot([1e-3,1e-1],[1e-13,1e-9])
            print(" J ext", self.ext_j)
            if do_all_plots:
                rlim = [1./2048, 0.1]
                axbonk(axd1[0][0], xlabel='r',ylabel=r'$\Omega$',xscale='log',yscale='log', xlim=rlim,ylim=self.ext_Omega.minmax)
                axbonk(axd1[1][0], xlabel='r',ylabel=r'$I_{ii}$',xscale='log',yscale='log',xlim=rlim, ylim=self.ext_I.minmax)
                axbonk(axd1[0][1], xlabel='r',ylabel=r'$j$',xscale='log',yscale='log',xlim=rlim, ylim=self.ext_j.minmax)
                axbonk(axd1[1][1], xlabel='r',ylabel=r'$\Omega/\rho$',xscale='log',yscale='log',xlim=rlim)
                outname = '%s/%s_vel_radial_to_look_at_core_c%04d'%(dl.output_directory,prefix,core_id)
                fig.savefig(outname)
                print("saved "+outname)



if 'do_all_plots' not in dir():
    do_all_plots = False
if 'clobber' not in dir():
    clobber=True


import three_loopers as tl

#if 'ang_tool1' not in dir():
#    ang_tool1=ang_tool(tl.looper1)
#    ang_tool1.run(do_all_plots=True,core_list=[51])


reload(trackage)
if 'looper_smol' not in dir() or clobber:
    looper_smol=looper.core_looper(directory=dl.sims['u05'])
    file_list = ['/data/cb1/Projects/P19_CoreSimulations/CoreSets/u05_every_ten/all_primitives_c0010.h5']
    for nfile,fname in enumerate(file_list):
        looper_smol.load_loop(fname)
ang_tool_smol = ang_tool(looper_smol)
ang_tool_smol.run()

if 1:
#if 'looper_smol' not in dir() or clobber:
#    looper_smol=looper.core_looper(directory=dl.sims['u05'])
#    file_list = ['/data/cb1/Projects/P19_CoreSimulations/CoreSets/u05_every_ten/all_primitives_c0010.h5']
#    for nfile,fname in enumerate(file_list):
#        looper_smol.load_loop(fname)
#ang_tool_smol = ang_tool(looper_smol)
#ang_tool_smol.run()
    if 'ang_tool1' not in dir() or clobber:
        ang_tool1=ang_tool(tl.looper1)
        ang_tool1.run(do_all_plots=do_all_plots)#,core_list=[51])
        ang_tool2=ang_tool(tl.looper2)
        ang_tool2.run(do_all_plots=do_all_plots)
        ang_tool3=ang_tool(tl.looper3)
        ang_tool3.run(do_all_plots=do_all_plots)
    for tool in [ang_tool1, ang_tool2, ang_tool3]:
        prefix = tool.this_looper.out_prefix
        plt.clf()
        ok = tool.times>0
        for core_id in tool.cores_used:
            plt.scatter( tool.times[ok], tool.slopes[core_id])
        plt.savefig('plots_to_sort/%s_Omega_slope.png'%(prefix))

    for tool in [ang_tool1, ang_tool2, ang_tool3]:
        prefix = tool.this_looper.out_prefix
        plt.clf()
        tool.times = tool.times[ tool.times > 0]
        these_slopes=[]
        slopes = np.zeros([tool.times.size, tool.cores_used.size])
        for nc, core_id in enumerate(tool.cores_used):
            slopes[:,nc] = tool.slopes[core_id]


        plt.clf()
        plt.boxplot(slopes.transpose())
        plt.savefig('plots_to_sort/%s_Omega_boxplot.pdf'%prefix)
        for core_id in tool.cores_used:
            plt.scatter( tool.times, tool.slopes[core_id])
        plt.savefig('plots_to_sort/%s_Omega_slope.png'%(prefix))
