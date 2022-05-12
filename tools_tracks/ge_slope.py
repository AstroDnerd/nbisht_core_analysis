from starter2 import *
import xtra_energy
import three_loopers_u500 as TL
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
sim_list=['u501']

def plot_phi(this_looper,core_list=None, do_plots=True):
    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    ge_array=[]
    gmm_array=[]
    frame = this_looper.target_frame
    for core_id in core_list:
        print('Potential %s %d'%(this_looper.sim_name,core_id))

        ds = this_looper.load(frame)
        xtra_energy.add_energies(ds)
        ms = trackage.mini_scrubber(this_looper.tr,core_id)
        c = nar([ms.mean_x[-1], ms.mean_y[-1],ms.mean_z[-1]])
        
        rsph = ds.arr(8.0/128,'code_length')
        sp = ds.sphere(c,rsph)

        GE = np.abs(sp['grav_energy'])
        dv = np.abs(sp['cell_volume'])
        RR =sp['radius']
        DD = sp['density']

        #2d distribution of GE vs r

        gbins = np.geomspace( GE[GE>0].min(), GE.max(),65)
        rbins = np.geomspace( RR [RR >0].min(), RR .max(),67)
        r_cen = 0.5*(rbins[1:]+rbins[:-1]) #we'll need this later.
        hist, xb, yb = np.histogram2d( RR , GE, bins=[rbins,gbins],weights=dv)

        #clear out annoying stragglers in the distribution.
        #any point that doesn't have any neighbors is eliminated.
        #h2 is the histogram to do math with.
        if 1:
            #h2 is the histogram.
            #we'll remove any stragglers.
            h2 = hist+0
            shifter = np.zeros(nar(h2.shape)+2)
            cuml = np.zeros(h2.shape)
            c_center = slice(1,-1)
            #hb is just "is h2 nonzero"
            #we'll slide this around to look for neighbors
            hb = (h2>0)
            shifter[c_center,c_center]=hb
            nx,ny=shifter.shape
            for i in [0,1,2]:
                for j in [0,1,2]:
                    if i==1 and j==1:
                        continue
                    s1 = slice(i,nx-2+i)
                    s2 = slice(j,ny-2+j)
                    cuml += shifter[s1,s2]
            #kill points that don't have neighbors.
            h2 *= (cuml >0)

        if 1:
            #Compute the upper bound of the histogram
            #smooth it
            #look for the point where the slope goes up.
            #but to avoid wiggles, it has to first come down a lot.

            #the upper bound of the distribution
            #compute upper_envelope
            y = np.arange(h2.shape[1])
            y2d = np.stack([y]*h2.shape[0])
            argmax = np.argmax(y2d*(h2>0),axis=1)
            upper_envelope = gbins[argmax]
            keepers=upper_envelope>1

            #smooth for wiggles.
            UE = gaussian_filter(upper_envelope[keepers],1)

            #
            # Find the inflection point where the slope comes up again.
            #
            #the slope
            DUE = UE[1:]-UE[:-1]
            #the max up to this radius
            cummax=np.maximum.accumulate( UE)
            #where the slope is positive
            ok = DUE > 0
            #and not too close to the center (for noise)
            ok = np.logical_and(ok , r_cen[keepers][1:]>1e-3)
            #and it has to come down below half its maximum
            ok = np.logical_and(ok, UE[1:]<0.5*cummax[1:])

            index = np.argmax(r_cen[keepers])
            if ok.any():
                #find out the radius where the inflection happens.
                index = np.where(ok)[0][0]
            R_KEEP = r_cen[keepers][index]

        #get the zones in the sphere that are
        #within R_KEEP
        ok_fit = np.logical_and(RR  < R_KEEP, GE>0)
        ok_fit = np.logical_and(ok_fit, RR  > r_cen[keepers].min())

        rok=RR[ok_fit].v

        #
        # Binding energy and GMM/R
        #
        ge_total = (GE[ok_fit]*dv[ok_fit]).sum()
        mtotal = (sp['cell_mass'][ok_fit]).sum()
        gmm      = G*mtotal**2/R_KEEP
        ge_array.append(ge_total)
        gmm_array.append(gmm)

        #Some plotting and fitting.
        if do_plots:
            if 1:
                #Just fit GE
                def plain_powerlaw( x, q, r0):
                    return q*x+r0
                popt, pcov=curve_fit(plain_powerlaw, np.log10(rok), np.log10(GE[ok_fit]))
                GE_fit_line=10**plain_powerlaw(np.log10(rok), *popt)

            if 1:
                #Fit density
                #Maybe its not necessary to histogram first, but it makes plotting easier.
                rbins = np.geomspace( RR [RR >0].min(), RR .max(),67)
                dbins = np.geomspace( DD[DD>0].min(), DD.max(),65)
                dhist, xbdr, ybdr = np.histogram2d( RR , DD, bins=[rbins,dbins],weights=dv)

                dok=DD[ok_fit]
                def powerlaw_r0_rkeep( r, q, rho0):
                    return q*np.log10(r/R_KEEP)+np.log10(rho0)
                poptd, pcovd=curve_fit(powerlaw_r0_rkeep, rok, np.log10(dok))

            fig,ax=plt.subplots(1,2)
            ax0=ax[0]; ax1=ax[1]

            if 1:
                #plot GE
                ax0.plot(r_cen[keepers], UE)
                ax0.scatter( r_cen[keepers][index],UE[index],c='r')
                pch.helper(h2,xb,yb,ax=ax0,transpose=False)
                ax0.plot( rok, GE_fit_line,c='r')
            if 1:
                #plot density
                pch.helper(dhist,xbdr,ybdr,ax=ax1,transpose=False)
                axbonk(ax1,xscale='log',yscale='log',xlabel='r',ylabel='rho')
                #density_fit_line=10**( poptd[0]*np.log10(rok)+np.log10(poptd[1]))
                density_fit_line=10**( powerlaw_r0_rkeep( rok, *poptd))
                ax1.plot( rok, density_fit_line,c='g')

                if 1:
                    rmin=r_cen[keepers].min()
                    ok2 = (RR < R_KEEP)*(RR >  rmin)
                    M = sp['cell_mass'][ok2].sum()
                    coe = 1/(8*np.pi*G)
                    power=2*poptd[0]+2
                    #print('DANS POWER',power)
                    #phi_del_squ_analy = (coe*G**2*M**2*rok**(power))/R_KEEP**(2*poptd[0]+6)
                    alpha=poptd[0]
                    rho0=poptd[1]
                    #phi_del_squ_analy = (4*np.pi*G*rho0*R_KEEP**(-alpha)*(2*alpha+5)/(alpha+3))**2*rok**power
                    phi_del_squ_analy = (4*np.pi*G*rho0*R_KEEP**(-alpha)*(alpha+2)/(alpha+3))**2*rok**power
                    ax0.plot( rok, phi_del_squ_analy ,c='g')

                    ratio=GE_fit_line/phi_del_squ_analy
                    #print(phi_del_squ_analy)
                    fig3,ax3=plt.subplots(1,1)
                    ax3.scatter(GE_fit_line,phi_del_squ_analy/GE_fit_line,c='k')
                    #axbonk(ax3,xscale='log',yscale='log',xlabel='GE',ylabel='d')
                    fig3.savefig('plots_to_sort/ratio.png')
                    #print('poot',(ratio).mean(), ratio.std())

                if 0:
                    #works pretty well
                    M = sp['cell_mass'].sum()
                    coe = 1/(8*np.pi*G)
                    power=2*poptd[0]+2
                    phi_del_squ_analy = (coe*G**2*M**2*rbins**(power))/RR.max()**(2*poptd[0]+6)
                    #print(phi_del_squ_analy)
                    ax0.plot( rbins, phi_del_squ_analy ,c='k')


            if 0:
                #color the upper envelope
                #to make sure we get it right.
                print(hist.shape)
                y = np.arange(hist.shape[1])
                y2d = np.stack([y]*hist.shape[0])
                argmax = np.argmax(y2d*(hist>0),axis=1)
                ind = np.arange( hist.shape[0])
                #take = np.ravel_multi_index(nar([argmax,ind]),hist.shape)
                take = np.ravel_multi_index(nar([ind,argmax]),hist.shape)
                h1=hist.flatten()
                h1[take]=hist.max()
                h1.shape=hist.shape
                pch.helper(h1,xb,yb,ax=ax0,transpose=False)

            outname='plots_to_sort/%s_c%04d_potfit'%(this_looper.sim_name,core_id)
            axbonk(ax0,xscale='log',yscale='log',xlabel='r',ylabel='grad phi sq')
            fig.savefig(outname)
            print(outname)
    return ge_array,gmm_array

if 'ge' not in dir() or True:
    for sim in TL.loops:
        if sim != 'u501':
            continue
        all_cores=np.unique( TL.loops[sim].tr.core_ids)
        core_list=list(all_cores)
        core_list=[323]
        ge,gmm=plot_phi( TL.loops[sim],core_list=core_list, do_plots=True)

fig,ax=plt.subplots(1,2)
ax0=ax[0];ax1=ax[1]
ax0.scatter(ge,gmm)
ax1.scatter(ge,nar(gmm)/nar(ge))
ax0.plot(ge,ge)
axbonk(ax0,xlabel='GE',ylabel='GMM/R',xscale='log',yscale='log')
axbonk(ax1,xlabel='GE',ylabel='GMM/R/GE', xscale='log')
fig.savefig('plots_to_sort/masses')
