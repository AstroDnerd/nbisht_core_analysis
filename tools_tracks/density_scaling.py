from starter2 import *
import data_locations as dl
import davetools
reload(davetools)
from scipy import signal
from scipy.optimize import curve_fit

plt.close('all')

def gauss_me(x, off,a, b):
    return off * np.exp(-(((x-a)**2)/(2*(b**2))))
def gauss_fit(centers,this_spectra, fit_fwhm_only=True):
    total = this_spectra.sum()
    vbar_est= (this_spectra*centers).sum()/total
    sigma2 = (this_spectra*(centers-vbar_est)**2).sum()/total
    sigma_est = (sigma2)**0.5
    norm_est = this_spectra.max()
    output = {'vbar_est':vbar_est,'total':total,'sigma_est':sigma_est, 'norm_est':norm_est}
    try:
        popt, pcov = curve_fit(gauss_me, centers, this_spectra)#, p0=[norm_est,vbar_est,sigma_est] )
        output.update( {'fit_norm':popt[0],'fit_center':popt[1],'fit_width':popt[2]})

    except:
        pass
    return output

extents = davetools.extents
file_list=glob.glob('%s/*h5'%dl.sixteen_frame)
file_list=file_list[:3]    

if 'this_looper' not in dir():
    this_looper=looper.core_looper(directory=dl.enzo_directory)
    for nfile,fname in enumerate(file_list):
        this_looper.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    thtr = this_looper.tr
    thtr.sort_time()
    all_cores = np.unique(thtr.core_ids)

core_list=all_cores
rm = rainbow_map(len(all_cores))

if 'rho_extents' not in dir():
    rho_extents=extents()
    x_extents=extents()
    r_extents=extents()
    alpha_extents=extents()
    for nc,core_id in enumerate(all_cores):
        ms = trackage.mini_scrubber(thtr,core_id)
        if ms.nparticles == 1:
            continue
        density = thtr.c([core_id],'density')
        t2 = np.tile(thtr.times,(ms.nparticles,1))
        alpha = density*t2**2
        r = ms.r
        x = r[t2>0]/t2[t2>0]
        rho_extents(density)
        alpha_extents(alpha[t2>0])
        r_extents(r)
        x_extents(x)

for nc,core_id in enumerate(core_list):
    asort =  np.argsort(thtr.times)
    if asort != sorted(asort):
        print("Warning: times not sorted.")
    n0=asort[0]
    tsorted = thtr.times[asort]
    fig=plt.figure(figsize=(8,8))
    axa=fig.subplots(2,2)
    #fig,axa=plt.subplots(1,2, figsize=(8,8))
    axd1=axa[0][0]
    axd2=axa[1][0]
    axd3=axa[0][1]
    axd4=axa[1][1]
    #fig2,ax666=plt.subplots(1,2)
    fig2=plt.figure(figsize=(8,8))
    ax666 = fig2.subplots(2,2)
    ax200 = ax666[0][0]
    ax201 = ax666[1][0]
    ax202 = ax666[0][1]
    ax203 = ax666[1][1]

    ms = trackage.mini_scrubber(thtr,core_id)
    if ms.nparticles == 1:
        continue
    density = thtr.c([core_id],'density')
    tmap=rainbow_map(ms.ntimes)
    norm = mpl.colors.Normalize()
    norm.autoscale( np.log10(density[:,n0]))
    cmap = mpl.cm.jet
    widths=[]
    means=[]
    timelist=[]

    alpha_mean_l=[]
    alpha_std=[]
    timelist_1=[]

    density_mean_l=[]
    density_std=[]

    r2_list=[]


    for n_count,n_time in enumerate(asort):
        time=thtr.times[n_time]
        if time == 0:
            continue
        c=tmap(n_count,ms.nparticles)
        this_r=ms.r[:,n_time]+0
        r_un = nar(sorted(np.unique(this_r)))

        this_x = this_r+0
        this_time = thtr.times[n_time]
        this_d =density[:,n_time] 
        alpha = this_d+0
        #if n_count == 0:
        #    alpha0 = alpha +0
        #    alpha_mean = alpha0.mean()

        #if this_time > 0:
        #   this_x *= this_time
        #   alpha *= (1-(this_time/0.046)**2)**2.1/alpha_mean
        #alpha = np.log(this_d)/np.log(this_time)
        alpha = np.log(this_d)/this_time**0.5 #flat then rises
        alpha = np.log(this_d)/this_time**0.5


        density_mean_l.append(np.mean(this_d))
        density_std.append(np.std(density))
        alpha_mean_l.append(np.mean(alpha))
        alpha_std.append(np.std(alpha))
        timelist_1.append(this_time)
        r2_list.append( np.sqrt( np.mean( this_r**2)))

        axd1.scatter(this_r,density[:,n_time],c=c,label=thtr.times[n_time],s=0.1)
        axd1.plot(r_un, 100*(r_un/1e-2)**-2,c='k',linewidth=0.1)
        axd2.scatter(this_x,alpha,c=c,label=thtr.times[n_time],s=0.1)
        axd3.scatter([this_time]*len(alpha),alpha,c=c,label=thtr.times[n_time],s=0.1)
        axd3.scatter([this_time+0.001]*len(this_d),this_d,c=[[0.3]*4]*len(this_d),label=thtr.times[n_time],s=0.1)

        #ax201.scatter(this_r,density[:,n_time],c=c,label=thtr.times[n_time],s=0.1)

        c=tmap(n_count)
        histout=ax201.hist( np.log10(this_d), color=c,histtype='step')
        hist, bins, wut=histout
        bc = 0.5*(bins[1:]+bins[:-1])
        fit_stuff=gauss_fit(bc, hist)
        #axd4.plot( bc, gauss_me( bc, fit_stuff['norm_est'], fit_stuff['vbar_est'], fit_stuff['sigma_est']))
        if 'fit_norm' in fit_stuff:
            axd4.plot( bc, gauss_me( bc, fit_stuff['fit_norm'], fit_stuff['fit_center'], fit_stuff['fit_width']),c=c)
            widths.append(fit_stuff['fit_width'])

            means.append(fit_stuff['fit_center'])
            timelist.append(this_time)

        x_un=nar(sorted(np.unique(this_x)))
        #ax201.plot(x_un, 1e-2*x_un**0,c='k',linewidth=0.1)
    davetools.axbonk(axd1,xscale='log',yscale='log',xlabel='r',ylabel=r'$\rho$',
                     xlim=r_extents.minmax, ylim=rho_extents.minmax)
    davetools.axbonk(axd2,xscale='log',yscale='log',xlabel='r',ylabel=r'$\rho$')
    #                 xlim=r_extents.minmax, ylim=rho_extents.minmax)
    davetools.axbonk(axd4,xscale='linear',yscale='log',xlabel='rho',ylabel='N',ylim=[1,50])
    davetools.axbonk(axd3,xscale='linear',yscale='log',xlabel='t',ylabel=r'$\alpha$')
                     #xlim=x_extents.minmax, ylim=alpha_extents.minmax)
    davetools.axbonk(ax201,xscale='linear',yscale='log',xlabel='rho',ylabel='N')
    outname = '/home/dcollins4096/PigPen/density_rescale_c%04d'%core_id
    fig.savefig(outname)
    print("saved "+outname)
    plt.close(fig)

    ax200.plot(timelist,widths,label='width')
    ax200.plot(timelist,means,label='means')
    ax202.plot(timelist_1,r2_list,marker='*')
    #ax203.plot(timelist_1,density_mean_l,marker='*',c='k')
    ax203.plot(timelist_1,alpha_mean_l,marker='*',c='r')
    ax203.set_yscale('log')
    import scatter_fit
    scatter_fit.scatter_fit(ax202,nar(timelist_1),nar(r2_list),log=False)
    ax200.legend(loc=0)
    #outname = 'image_tracks/density_scale_c%04d'%core_id
    outname = '/home/dcollins4096/PigPen/hist_prop_c%04d'%core_id
    fig2.savefig(outname)
    print("saved",outname)
    plt.close(fig2)
