from starter2 import *
import data_locations as dl
import davetools
reload(davetools)
import scatter_fit
reload(scatter_fit)

plt.close('all')


file_list=glob.glob('%s/*h5'%dl.sixteen_frame)
#for debug purposes you may want a reduced list 
#file_list=file_list[:3]    

if 'this_looper' not in dir():
    this_looper=looper.core_looper(directory=dl.enzo_directory)
    for nfile,fname in enumerate(file_list):
        this_looper.load_loop(fname)
        print( "Reading file %d of %d"%(nfile,len(file_list)))
    thtr = this_looper.tr
    thtr.sort_time()
    all_cores = np.unique(thtr.core_ids)

core_list=all_cores
rm = rainbow_map(len(all_cores))

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

r_extents.minmax[0]=0.5/2048
asort =  np.argsort(thtr.times)
if (asort != sorted(asort)).any():
    print("Warning: times not sorted.")
do_all_plots=False

from scipy.optimize import curve_fit
def powerlaw(r,rho0, r0, alpha):
    return alpha*np.log10(r/r0) + np.log10(rho0)

if 'alpha' not in dir() or True:
    alpha = []
    mass = []
    density_0=[]
    mean_div = []
    alpha_dict={}
    ft_lst_rho0=[]
    ft_lst_r0=[]
    ft_lst_alpha=[]
    mini_alphas={}


    for nc,core_id in enumerate(core_list):

        #miniscrubber computes distance, r^2, several other quantities
        ms = trackage.mini_scrubber(thtr,core_id)
        if ms.nparticles == 1:
            continue

        density = thtr.c([core_id],'density')


        #Compute alpha
        #compute mean initial density.
        this_r = ms.r.flatten()
        min_r = 1./2048
        this_r[ this_r < min_r ] = min_r    

        this_rho = density.flatten()
        popt, pcov = curve_fit(powerlaw, this_r, np.log10(this_rho), p0=[1,1,-2])
        fit_rho0, fit_r0, fit_alpha = popt
        ft_lst_rho0.append(fit_rho0)
        ft_lst_r0.append(fit_r0)
        ft_lst_alpha.append(fit_alpha)

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
            mini_alphas[core_id]=[]
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
                mini_alphas[core_id].append(lilfit[0])


            title=""
            title+=r'$\alpha = %0.3f\ \rho_0 = $%s$ r_0=$%s'%(alpha[-1], expform(fit_rho0),expform(fit_r0))
            axd1.set_title(title)
            axd2.scatter(time_list_dumb, mini_alphas[core_id],c=color_list_dumb)
            axd2.plot(time_list_dumb, mini_alphas[core_id],c=[0.5]*3)
            axd2.plot(time_list_dumb, [fit_alpha]*len(time_list_dumb),c='k')

            axbonk(axd1,xscale='log',yscale='log',xlabel='r',ylabel=r'$\rho$',
                             xlim=r_extents.minmax, ylim=rho_extents.minmax)
            axbonk(axd2,xlabel=r'$t$',ylabel=r'$\alpha(t)$',ylim=[-4,2])
            #axbonk(axd2,xscale='log',yscale='log',xlabel='r',ylabel=r'$\rho$',
            #                 xlim=r_extents.minmax, ylim=rho_extents.minmax)
            axd1.set_ylim(rho_extents.minmax)
            axd1.set_xlim(r_extents.minmax)
            print("R EXTENTS", r_extents.minmax)
            outname = '%s/density_radius_c%04d'%(dl.output_directory,core_id)
            fig.savefig(outname)
            print("saved "+outname)
            plt.close(fig)

if 0:
    mass=nar(mass)
    fig, axes=plt.subplots(1,2)
    ax0 = axes[0]; ax1=axes[1]
    ax0.scatter(density_0,alpha)
    ax1.scatter(mass,alpha)
    davetools.axbonk(ax0,ylabel=r'$\alpha$', xlabel=r'$\langle \rho \rangle$', xscale='log',yscale='linear')
    davetools.axbonk(ax1,ylabel=r'$\alpha$', xlabel=r'$M(t=0)$', xscale='log',yscale='linear')
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


"""
if 'alpha_proxy' not in dir() or True:
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


"""
