from starter2 import *
from scipy.optimize import curve_fit
import data_locations as dl
import davetools
reload(davetools)
import scatter_fit
reload(scatter_fit)

plt.close('all')

from three_loopers import *

this_simname = 'u05'
#for debug purposes you may want a reduced list 
#file_list=file_list[:3]    

def powerlaw(r,rho0, r0, alpha):
    return alpha*np.log10(r/r0) + np.log10(rho0)

class alpha_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper

        self.rho_extents=davetools.extents()
        self.r_extents=davetools.extents()
        self.alpha = []
        self.mass = []
        self.density_0=[]
        self.mean_div = []
        self.alpha_dict={}
        self.ft_lst_rho0=[]
        self.ft_lst_r0=[]
        self.ft_lst_alpha=[]
        self.mini_alphas={}
        self.collapse_times=[]
        self.cores_used=[]

    def run(self, core_list=None, do_all_plots=False):
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list=all_cores
        rm = rainbow_map(len(all_cores))

        if len(self.rho_extents.minmax) == 0:
            for nc,core_id in enumerate(all_cores):
                ms = trackage.mini_scrubber(thtr,core_id)
                if ms.nparticles == 1:
                    continue
                density = thtr.c([core_id],'density')
                self.rho_extents(density)
                self.r_extents(ms.r)

        self.r_extents.minmax[0]=0.5/2048
        asort =  np.argsort(thtr.times)
        if (asort != sorted(asort)).any():
            print("Warning: times not sorted.")

        for nc,core_id in enumerate(core_list):

            #miniscrubber computes distance, r^2, several other quantities
            ms = trackage.mini_scrubber(thtr,core_id)
            if ms.nparticles < 3:
                continue
            self.cores_used.append(core_id)


            density = thtr.c([core_id],'density')

            if 1:
                #collapse time
                if (density>1e3).any():
                    first_collapse_index = np.where( (density>1e3).sum(axis=0) >0 )[0][0] 
                    t_collapse = thtr.times[first_collapse_index]
                    self.collapse_times.append(t_collapse)
                else:
                    t_collapse = -1
                    self.collapse_times.append(t_collapse)

            #Compute alpha
            #compute mean initial density.
            this_r = ms.r.flatten()
            min_r = 1./2048
            this_r[ this_r < min_r ] = min_r    
            fit = scatter_fit.scatter_fit(plt,this_r,density.flatten())
            self.alpha.append(fit['fit'][0])
            self.alpha_dict[core_id]=self.alpha[-1]
            self.density_0.append( density[:,0].mean()) #this is not good

            this_rho = density.flatten()
            popt, pcov = curve_fit(powerlaw, this_r, np.log10(this_rho), p0=[1,1,-2])
            fit_rho0, fit_r0, fit_alpha = popt
            self.ft_lst_rho0.append(fit_rho0)
            self.ft_lst_r0.append(fit_r0)
            self.ft_lst_alpha.append(fit_alpha)

            #compute initial mass.
            for nf in [0]:
                dx=1./2048
                nx=2048
                x =np.floor(thtr.c([core_id],'x')/dx)[:,nf]
                y =np.floor(thtr.c([core_id],'y')/dx)[:,nf]
                z =np.floor(thtr.c([core_id],'z')/dx)[:,nf]
                this_density = density[:,nf]
                cell_volume = thtr.c([core_id],'cell_volume')[:,nf]
                index = x + nx*(y * nx*z)
                ar = np.argsort(index)
                rs = np.argsort(ar)
                isorted=index[ar]
                mask = np.ones_like(this_density,dtype='bool')
                mask[1:] = isorted[1:]-isorted[:-1] != 0
                mask2 = mask[ rs]
                self.mass.append((this_density[mask2]*cell_volume[mask2]).sum())

            if do_all_plots:
                fig, axes=plt.subplots(1,2)
                axd1 = axes[0]; axd2 = axes[1]
                #other_fit = other_fit(this_r, density.flatten())
                this_r = np.array([min_r, 0.3])
                #this_r[ this_r < min_r] = min_r
                axd1.plot( this_r, 10**powerlaw(this_r, popt[0],popt[1],popt[2]))
                #axd1.plot( this_r, 0.1*(this_r/0.3)**alpha[-1],c='k')
                #axd1.plot( this_r, 0.1*(this_r/0.3)**-0.5,c=[0.9]*3)
                #axd1.plot( this_r, 0.1*(this_r/0.3)**-1.0,c=[0.9]*3)
                #axd1.plot( this_r, 0.1*(this_r/0.3)**-1.5,c=[0.9]*3)
                #axd1.plot( this_r, 0.1*(this_r/0.3)**-2.0,c=[0.9]*3)
                #axd1.plot( this_r, 0.1*(this_r/0.3)**-2.5,c=[0.9]*3)
                tmap=rainbow_map(ms.ntimes)
                self.mini_alphas[core_id]=[]
                time_list_dumb=[]
                color_list_dumb=[]

                for n_count,n_time in enumerate(asort):
                    time=thtr.times[n_time]
                    if time == 0:
                        continue
                    c=tmap(n_count,ms.nparticles)
                    time_list_dumb.append(time)
                    color_list_dumb.append(c[0])
                    this_r=ms.r[:,n_time]+0
                    this_r[ this_r < min_r] = min_r
                    this_density=density[:,n_time]
                    #popt, pcov = curve_fit(powerlaw, this_r, np.log10(this_density), p0=[1,1,-2])
                    lilfit=np.polyfit( np.log10(this_r), np.log10(this_density), 1)
                    axd1.scatter(this_r,this_density,c=c,label=thtr.times[n_time],s=0.1)
                    axd1.plot( this_r, 10**(lilfit[0]*np.log10(this_r)+lilfit[1]),c=c[0])
                    self.mini_alphas[core_id].append(lilfit[0])


                title=""
                title+=r'$\alpha = %0.3f\ \rho_0 = $%s$ r_0=$%s'%(self.alpha[-1], expform(fit_rho0),expform(fit_r0))
                axd1.set_title(title)
                axbonk(axd1,xscale='log',yscale='log',xlabel='r',ylabel=r'$\rho$',
                                 xlim=self.r_extents.minmax, ylim=self.rho_extents.minmax)
                axd1.set_ylim(self.rho_extents.minmax)
                axd1.set_xlim(self.r_extents.minmax)

                axd2.scatter(time_list_dumb, self.mini_alphas[core_id],c=color_list_dumb)
                axd2.plot(time_list_dumb, self.mini_alphas[core_id],c=[0.5]*3)
                #axd2.plot(time_list_dumb, [fit_alpha]*len(time_list_dumb),c='k')
                axd2.plot(time_list_dumb, [-1]*len(time_list_dumb),  c=[0.5]*4)
                axd2.plot(time_list_dumb, [-0.5]*len(time_list_dumb),c=[0.5]*4)
                axd2.plot(time_list_dumb, [-1.5]*len(time_list_dumb),c=[0.5]*4)
                axd2.plot(time_list_dumb, [-2]*len(time_list_dumb),  c=[0.5]*4)

                axbonk(axd2,xlabel=r'$t$',ylabel=r'$\alpha(t)$',ylim=[-4,2])
                #axbonk(axd2,xscale='log',yscale='log',xlabel='r',ylabel=r'$\rho$',
                #                 xlim=r_extents.minmax, ylim=rho_extents.minmax)
                outname = '%s/%s_density_radius_c%04d'%(dl.output_directory,self.this_looper.out_prefix,core_id)
                fig.savefig(outname)
                print("saved "+outname)
                plt.close(fig)


    def do_plots(self):
            mass=nar(self.mass)
            sim_color={'u05':'r','u10':'g','u11':'b'}[self.this_looper.out_prefix]
            fig, axes=plt.subplots(2,2)
            ax0 = axes[0][0]; ax1=axes[0][1]
            ax2 = axes[1][0]; ax3=axes[1][1]
            ax0.scatter(self.density_0,self.alpha,c=sim_color)
            ax1.scatter(mass,self.alpha,c=sim_color)

            fig2, axes2=plt.subplots(1,1)
            bx0 = axes2
            if 'collapse_times' in dir(): #collapse_time can be found from density_time.py
                ok = nar(collapse_times) > 0
                t_ok=nar(collapse_times)[ok]
                ax2.scatter(t_ok, nar(self.alpha)[ok])
                ax3.scatter(t_ok, nar(self.mass)[ok],c=nar(sim_color)[ok])
                bx0.scatter(t_ok, nar(self.density_0)[ok])
                #ax3.scatter(t_ok, nar(density_0)[ok],c=nar(sim_color)[ok])
                pdb.set_trace()

            davetools.axbonk(ax0,xlabel=r'$\langle \rho \rangle$',ylabel=r'$\alpha$',  xscale='log',yscale='linear')
            davetools.axbonk(ax1,xlabel=r'$M(t=0)$',              ylabel=r'$\alpha$',   xscale='log',yscale='linear')
            davetools.axbonk(ax2,xlabel=r'$t_c$',                 ylabel=r'$\alpha$', xscale='log',yscale='linear')
            davetools.axbonk(ax3,xlabel=r'$t_c$',                 ylabel=r'$M(t=0)$',  xscale='log',yscale='log')
            davetools.axbonk(bx0,xlabel=r'$t_c$',   ylabel=r'$\langle \rho \rangle$',  xscale='log',yscale='linear')
            ax3.set_xlim([t_ok.min(),t_ok.max()])
            ax3.set_ylim(self.mass.min(),self.mass.max())
            ax1.set_xlim(self.mass.min(),self.mass.max())
            outname = '%s/alpha_mass.png'%dl.output_directory
            fig.savefig(outname)
            outname2 = '%s/tc_rho.png'%dl.output_directory
            fig2.savefig(outname2)
            print(outname)

    if 0:
        mass=nar(mass)
        fig, axes=plt.subplots(2,2)
        ax0 = axes[0][0]; ax1=axes[0][1]
        ax2 = axes[1][0]; ax3=axes[1][1]
        ax0.scatter(density_0,alpha)
        ax1.scatter(mass,alpha)
        #ax2.scatter(total_volume,alpha)
        #ax3.scatter(total_volume,mass)
        davetools.axbonk(ax0,ylabel=r'$\alpha$', xlabel=r'$\langle \rho \rangle$', xscale='log',yscale='linear')
        davetools.axbonk(ax1,ylabel=r'$\alpha$', xlabel=r'$M(t=0)$', xscale='log',yscale='linear')
        #vol_lim=[min(total_volume),max(total_volume)]
        mass_lim=[min(mass),max(mass)]
        davetools.axbonk(ax2,ylabel=r'$\alpha$', xlabel=r'$V(t=0)$', xscale='log',yscale='linear',xlim=vol_lim)
        davetools.axbonk(ax3,ylabel=r'$M(t=0)$', xlabel=r'$V(t=0)$', xscale='log',yscale='log',xlim=vol_lim,ylim=mass_lim)
        ax1.set_xlim(mass.min(),mass.max())
        outname = '%s/alpha_mass.png'%dl.output_directory
        fig.savefig(outname)
        print(outname)


#plt.clf()
#alpha=np.array(alpha)
#ft_lst_alpha=np.array(ft_lst_alpha)
#plt.scatter(alpha,1-(alpha/ft_lst_alpha))
#plt.savefig('plots_to_sort/alpha_alpha.png')
#plt.close('all')

#fig,ax=plt.subplots(1,1)
#ax.scatter(ft_lst_rho0,ft_lst_r0)
#axbonk(ax,xlabel=r'$\rho_0$',ylabel=r'$r_0$',xscale='log',yscale='log')
#fig.savefig('plots_to_sort/rho0_r0.png')


    if 'alpha_proxy' not in dir() and False:
        do_all_plots=True
        div_v=[]
        v_over_r=[]
        divv_extents=davetools.extents()
        divv_extents( thtr.track_dict['velocity_divergence'])
        proxy_extents=davetools.extents()
        for nc,core_id in enumerate(core_list):
            #miniscrubber computes distance, r^2, several other quantities
            ms = trackage.mini_scrubber(thtr,core_id)
            if ms.nparticles == 1:
                continue
            density = thtr.c([core_id],'density')
            div_v = thtr.c([core_id],'velocity_divergence')
            mag_v = thtr.c([core_id],'velocity_magnitude')
            divv_extents(div_v)
            if do_all_plots:
                tmap=rainbow_map(ms.ntimes)
                fig, axes=plt.subplots(2,2)
                ax00=axes[0][0]; ax01=axes[0][1]
                ax10=axes[1][0]; ax11=axes[1][1]

                for n_count,n_time in enumerate(asort):
                    time=thtr.times[n_time]
                    if time == 0:
                        continue
                    c=tmap(n_count,ms.nparticles)
                    this_r=ms.r[:,n_time]+0

                    ax00.scatter(this_r,div_v[:,n_time],c=c,label=thtr.times[n_time],s=0.1)
                    ax01.scatter(this_r,mag_v[:,n_time],c=c,label=thtr.times[n_time],s=0.1)


                    alpha_proxy = div_v[:,n_time]/mag_v[:,n_time]*ms.r[:,n_time]
                    proxy_extents(alpha_proxy)
                    ax10.scatter(this_r,alpha_proxy.flatten(),c=c,label=thtr.times[n_time],s=0.1)
                davetools.axbonk(ax00,xscale='log',yscale='linear',xlabel='r',ylabel=r'$\nabla\cdot v$',
                                 xlim=r_extents.minmax, ylim=rho_extents.minmax)
                ax00.set_yscale('symlog',linthreshy=1)
                ax00.set_ylim(-5e5,5e5)

                davetools.axbonk(ax01,xscale='log',yscale='log',xlabel='r',ylabel=r'$|v|$',
                                 xlim=r_extents.minmax, ylim=rho_extents.minmax)
                ax01.set_yscale('symlog',linthreshy=1)
                ax01.set_ylim(0,100)

                ax10.plot([1e-3,1e-1], [alpha_dict[core_id],alpha_dict[core_id]],c='k')
                ax10.plot([1e-3,1e-1], [-0.5,-0.5],c=[0.5]*4)
                ax10.plot([1e-3,1e-1], [-3,-3],c=[0.5]*4)
                davetools.axbonk(ax10,xscale='log',yscale='linear',xlabel='r',ylabel=r'$\nabla\cdot v/(v/r)$',
                                 xlim=r_extents.minmax, ylim=rho_extents.minmax)
                ax10.set_yscale('symlog',linthreshy=0.1)
                ax10.set_ylim(proxy_extents.minmax)

                outname = '%s/divergence_proxy_c%04d'%(dl.output_directory,core_id)
                fig.savefig(outname)
                print("saved "+outname)
                plt.close(fig)


import three_loopers as tl
do_all_plots=False
if 'clobber' not in dir():
    clobber=False
if 't1' not in dir() or clobber:
    t1 =alpha_tool(tl.looper1)
    t2 =alpha_tool(tl.looper2)
    t3 =alpha_tool(tl.looper3)
    t1.run(do_all_plots=do_all_plots)
    t2.run(do_all_plots=do_all_plots)
    t3.run(do_all_plots=do_all_plots)

import tools_tracks.density_time as DT
reload(DT)
if 'density_tool1' not in dir():
    density_tool1=DT.trial(tl.looper1)
    density_tool1.run(do_all_plots=False)
    density_tool2=DT.trial(tl.looper2)
    density_tool2.run(do_all_plots=False)
    density_tool3=DT.trial(tl.looper3)
    density_tool3.run(do_all_plots=False)

if 1:
    sim_color={'u05':'r','u10':'g','u11':'b'}#[this_tool.this_looper.out_prefix]
    fig, axes=plt.subplots(2,2)
    ax0 = axes[0][0]; ax1=axes[0][1]
    ax2 = axes[1][0]; ax3=axes[1][1]


    fig2, axes2=plt.subplots(1,1)
    bx0 = axes2

    density_tools = [density_tool1, density_tool2, density_tool3]
    alpha_tools = [t1, t2, t3]
    
    mass_ext = extents()
    for this_dt, this_at in zip(density_tools,alpha_tools):
        if 1:
            ax0.clear()
            ax1.clear()
            ax2.clear()
            ax3.clear()

        ax0.scatter(this_at.density_0,this_at.alpha,c=sim_color[this_at.this_looper.out_prefix])
        ax1.scatter(this_at.mass,this_at.alpha,c=sim_color[this_at.this_looper.out_prefix])
        print('col',len(this_dt.collapse_times))
        print('alph', len(this_at.alpha))
        c=sim_color[this_at.this_looper.out_prefix]
        ok = nar(this_dt.collapse_times) > 0
        t_ok=nar(this_dt.collapse_times)[ok]
        ax2.scatter(t_ok, nar(this_at.alpha)[ok],c=c)
        ax3.scatter(t_ok, nar(this_at.mass)[ok],c=c)
        bx0.scatter(t_ok, nar(this_at.density_0)[ok])
        #ax3.scatter(t_ok, nar(density_0)[ok],c=nar(sim_color)[ok])
        mass_ext(nar(this_at.mass))

        davetools.axbonk(ax0,xlabel=r'$\langle \rho \rangle$',ylabel=r'$\alpha$',  xscale='log',yscale='linear')
        davetools.axbonk(ax1,xlabel=r'$M(t=0)$',              ylabel=r'$\alpha$',   xscale='log',yscale='linear')
        davetools.axbonk(ax2,xlabel=r'$t_c$',                 ylabel=r'$\alpha$', xscale='log',yscale='linear')
        davetools.axbonk(ax3,xlabel=r'$t_c$',                 ylabel=r'$M(t=0)$',  xscale='log',yscale='log')
        davetools.axbonk(bx0,xlabel=r'$t_c$',   ylabel=r'$\langle \rho \rangle$',  xscale='log',yscale='linear')
        ax3.set_xlim([t_ok.min(),t_ok.max()])
        ax3.set_ylim(mass_ext.minmax)
        ax1.set_xlim(mass_ext.minmax)
        outname = '%s/alpha_mass.png'%dl.output_directory
        outname = '%s/%s_alpha_mass.png'%(dl.output_directory,this_dt.this_looper.out_prefix)
        fig.savefig(outname)
        outname2 = '%s/tc_rho.png'%dl.output_directory
        fig2.savefig(outname2)
        print(outname)
