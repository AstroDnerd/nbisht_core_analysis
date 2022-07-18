from starter2 import *
import xtra_energy
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
sim_list=['u501']

from collections import defaultdict
import r_inflection
#import three_loopers_u500 as TL
import three_loopers_six as TL
if 'inflection' not in dir():
    inflection = {}
for sim in TL.loops:
    if sim not in inflection:
        inflection[sim]=r_inflection.R_INFLECTION( TL.loops[sim])
        inflection[sim].run()

class slope_tool():
    def __init__(self, this_looper, inflection):
        self.this_looper=this_looper
        self.output = defaultdict(list)
        self.rinflection=inflection.rinflection
        self.rinflection_list=inflection.rinflection_list

    def run(self,core_list=None, do_plots=True,do_proj=True):
        print("CORE LIST",core_list)
        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        this_looper=self.this_looper
        frame = this_looper.target_frame
        ds = this_looper.load(frame)
        G = ds['GravitationalConstant']/(4*np.pi)
        xtra_energy.add_energies(ds)
        for core_id in core_list:
            print('GE slope %s %d'%(this_looper.sim_name,core_id))

            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            c = nar([ms.mean_x[-1], ms.mean_y[-1],ms.mean_z[-1]])
            
            R_SPHERE = 8/128
            rsph = ds.arr(R_SPHERE,'code_length')
            sp = ds.sphere(c,rsph)




            GE = np.abs(sp['grav_energy'])
            dv = np.abs(sp['cell_volume'])
            RR = sp['radius']
            DD = sp['density']

            R_KEEP = self.rinflection[core_id]

            #2d distribution of GE vs r

            gbins = np.geomspace( GE[GE>0].min(), GE.max(),65)
            rbins = np.geomspace( RR [RR >0].min(), RR .max(),67)
            r_cen = 0.5*(rbins[1:]+rbins[:-1]) #we'll need this later.
            hist, xb, yb = np.histogram2d( RR , GE, bins=[rbins,gbins],weights=dv)


            #get the zones in the sphere that are
            #within R_KEEP
            ok_fit = np.logical_and(RR  < R_KEEP, GE>0)

            rok=RR[ok_fit].v

            #
            # Binding energy and GMM/R
            #
            ge_total = (GE[ok_fit]*dv[ok_fit]).sum()
            mtotal = (sp['cell_mass'][ok_fit]).sum()
            gmm      = G*mtotal**2/R_KEEP
            self.output['ge_total'].append(ge_total)
            self.output['gmm'].append(gmm)

            if 1:
                #store stuff
                self.output['mass'].append(mtotal)
                self.output['r0'].append(R_KEEP)

            if 0:
                #Just fit GE
                def plain_powerlaw( x, q, r0):
                    return q*x+r0
                popt, pcov=curve_fit(plain_powerlaw, np.log10(rok), np.log10(GE[ok_fit]))
                GE_fit_line=10**plain_powerlaw(np.log10(rok), *popt)
                ouput['alpha_ge'].append(popt[0])
            if 1:
                #Better ge fit
                def plain_powerlaw( x, q, norm):
                    return (2*q+2)*x+norm#np.log10(G*mtotal/R_KEEP)
                popt, pcov=curve_fit(plain_powerlaw, np.log10(rok/R_KEEP), np.log10(GE[ok_fit]))
                A = 10**popt[1]
                alpha = popt[0]
                rmax = (rok).max()
                rmin = (rok).min()

                rfit = np.linspace(rmin,rmax,128)
                rr = 0.5*(rfit[1:]+rfit[:-1])
                dr = rfit[1:]-rfit[:-1]

                GE_fit_line=10**plain_powerlaw(np.log10(rr/R_KEEP), *popt)
                #GE_fit_line=10**plain_powerlaw(np.log10(rok/R_KEEP), *popt)
                self.output['alpha_ge'].append(alpha)
                self.output['A'].append(A)
                self.output['AA'].append( G**2*mtotal**2/R_KEEP**(2*alpha+2))

                R0=R_KEEP
                self.output['analytic'].append( 4*np.pi*A/( (R0**(2*alpha+2)*(2*alpha+5)))*(rmax**(2*alpha+5)-rmin**(2*alpha+5)))
                self.output['rmin'].append(rmin)
                self.output['rmax'].append(rmax)
                self.output['analytic2'].append( 4*np.pi*A/( (R0**(2*alpha+2)*(2*alpha+5)))*(R_KEEP**(2*alpha+5)))
                R0 = R_KEEP
                M = mtotal 
                E1 = (-G*M*M/R0**(2*alpha+6)/(2*alpha+5))*(R0**(2*alpha+5)-rmin**(2*alpha+5))

                self.output['ann_good'].append( E1)

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
                self.output['alpha_rho'].append(poptd[0])

            if do_proj:
                proj_axis=0
                dx = 1/2048
                nx = 2*R_SPHERE/dx

                pd = ds.proj('density',proj_axis,data_source=sp, center=c)
                frb=pd.to_frb(2*R_SPHERE,nx,center=c)
                fig2,rx=plt.subplots(1,2, figsize=(12,8))
                rx0=rx[0]; rx1=rx[1]
                column_density=frb['density']
                norm = mpl.colors.LogNorm( column_density.min(),column_density.max())
                rx0.imshow(column_density,norm=norm)

                column_density=frb[YT_grav_energy]
                linthresh=1
                norm = mpl.colors.SymLogNorm(linthresh, vmin=column_density.min(),vmax=0)
                rx1.imshow(column_density,norm=norm)

                center = nar(column_density.shape)/2
                circle = plt.Circle( center, R_KEEP/dx,fill=False)
                rx0.add_artist(circle)
                circle = plt.Circle( center, R_KEEP/dx,fill=False)
                rx1.add_artist(circle)
                fig2.savefig('plots_to_sort/cPhiProj_%s_c%04d'%(this_looper.sim_name, core_id))
            #Some plotting and fitting.
            if do_plots:

                fig,ax=plt.subplots(1,2)
                ax0=ax[0]; ax1=ax[1]

                if 0:
                    #plot GE
                    ax0.plot(r_cen[keepers], UE)
                    pch.helper(h2,xb,yb,ax=ax0,transpose=False)
                    #ax0.plot( rok, GE_fit_line,c='r')
                    ax0.plot( rr, GE_fit_line,c='r')
                    ax0.scatter( R_KEEP,UE[index],c='r')
                    AlsoA=colors.G*mtotal**2/(4*np.pi)*R_KEEP**(-4)
                    ax0.scatter( R_KEEP, A, c='g',marker='*')
                    ax0.scatter( R_KEEP, AlsoA, c='b',marker='*')
                if 1:
                    #plot density
                    pch.helper(dhist,xbdr,ybdr,ax=ax1,transpose=False)
                    axbonk(ax1,xscale='log',yscale='log',xlabel='r',ylabel='rho')
                    #density_fit_line=10**( poptd[0]*np.log10(rok)+np.log10(poptd[1]))
                    density_fit_line=10**( powerlaw_r0_rkeep( rok, *poptd))
                    ax1.plot( rok, density_fit_line,c='g')

                    if 0:
                        rmin=r_cen[keepers].min()
                        ok2 = (RR < R_KEEP)*(RR >  rmin)
                        M = sp['cell_mass'][ok2].sum()
                        coe = 1/(8*np.pi*G)
                        power=2*poptd[0]+2
                        #print('DANS POWER',power)
                        #phi_del_squ_analy = (coe*G**2*M**2*rok**(power))/R_KEEP**(2*poptd[0]+6)
                        #horse around
                        DanA = (coe*G**2*M**2)/R_KEEP**(2*poptd[0]+6)
                        phi_del_squ_analy = DanA*rok**(power)
                        self.output['DanA'].append( DanA)
                        alpha=poptd[0]
                        rho0=poptd[1]
                        #phi_del_squ_analy = (4*np.pi*G*rho0*R_KEEP**(-alpha)*(2*alpha+5)/(alpha+3))**2*rok**power
                        #phi_del_squ_analy = (4*np.pi*G*rho0*R_KEEP**(-alpha)*(alpha+2)/(alpha+3))**2*rok**power
                        ax0.plot( rok, phi_del_squ_analy ,c='g')

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


if 'stuff' not in dir():# or True:
    stuff={}
    for sim in ['u603']:
        #stuff={}
        all_cores=np.unique( TL.loops[sim].tr.core_ids)
        #core_list=list(all_cores)
        core_list=all_cores
        #core_list=[323]
        #stuff[sim]=plot_phi( TL.loops[sim],core_list=core_list, do_plots=False)
        stuff[sim]=slope_tool(TL.loops[sim], inflection[sim])
        stuff[sim].run(core_list=core_list,do_plots=False, do_proj=False)

if 1:
    for sim in stuff:
        ge   =nar(stuff[sim].output['ge_total'])
        gmm1 =nar(stuff[sim].output['gmm'])
        bgd  =nar(stuff[sim].output['alpha_ge'])
        #ann  =nar(stuff[sim]['analytic'])
        ann1=np.abs(nar(stuff[sim].output['analytic']))
        ann2 =np.abs(nar(stuff[sim].output['ann_good']))
        A    =nar(stuff[sim].output['A'])
        mass =nar(stuff[sim].output['mass'])
        r0   =nar(stuff[sim].output['r0'])
        rmin   =nar(stuff[sim].output['rmin'])
        rmax   =nar(stuff[sim].output['rmax'])
        alpha=nar(stuff[sim].output['alpha_ge'])
        alpha_d=nar(stuff[sim].output['alpha_rho'])
        AA1   =nar(stuff[sim].output['AA'])
        AA2 = (colors.G*mass**2/r0**(2*alpha+6))/(2*alpha+5)
        GG2 = AA2*(r0**(2*alpha+5)-rmin**(2*alpha+5))
        GG3 = 4*np.pi*A*(rmax**(2*alpha+5)-rmin**(2*alpha+5))/(r0**(2*alpha+5)*(2*alpha+5))
        AA3 =  nar(stuff[sim].output['DanA'])

        agd = bgd #(nar(bgd)-2)/2

        gmm2 = gmm1/(4*agd+10)
        #gmm2 = gmm1*(3+agd)/(5+2*agd)

        gmm=gmm1
        fig,ax=plt.subplots(4,2)
        ax0=ax[0][0];ax1=ax[0][1]
        ax2=ax[1][0];ax3=ax[1][1]
        ax4=ax[2][0];ax5=ax[2][1]
        ax6=ax[3][0];ax7=ax[3][1]

        ax0.scatter(ge,gmm)
        ax0.plot(ge,ge)
        axbonk(ax0,xlabel='GE',ylabel='GMM/R',xscale='log',yscale='log')


        ax1.scatter(ge,nar(ge)/nar(gmm))
        axbonk(ax1,xlabel='GE',ylabel='GE/GMM/R', xscale='log')


        ax2.scatter( ge, ann1)
        ax2.plot(ge,ge)
        axbonk(ax2,xlabel='GE',ylabel=r'$\phi=A(r/r0)^(2\alpha+2)$',xscale='log',yscale='log')

        ann3=4*np.pi*A/(r0**(2*alpha+2)*(2*alpha+5))*(r0**(2*alpha+5)-rmin**(2*alpha+5))
        #ax3.scatter(4*np.pi*A,AA2, c='g')
        AlsoA=colors.G*mass**2/(4*np.pi)*r0**(-4)
        ax3.scatter(A, AlsoA)
        ax3.plot(A,A)

        #ax3.plot( 4*np.pi*A, 4*np.pi*A/64)

        axbonk(ax3,xscale='log',yscale='log',xlabel='A',ylabel='GM^2/(4 pi r^4)')


        ax4.scatter(ge,ann2)
        ax4.plot(ge,ge)
        axbonk(ax4,xlabel='GE',ylabel='M^r^(2alpha_phi+2)',xscale='log',yscale='log')

        ax5.scatter( alpha,alpha_d)
        alpha_ext=extents(alpha)
        alpha_ext(alpha_d)
        axbonk(ax5,xlabel='alpha-phi',ylabel='alpha-d',xlim=alpha_ext.minmax,ylim=alpha_ext.minmax)

        
        ann5  = (colors.G*mass**2/r0**(2*alpha_d+6)/(2*alpha_d+5))*(r0**(2*alpha_d+5)-rmin**(2*alpha_d+5))
        #ann5b  = (colors.G*mass**2/r0**(2*alpha_d+6)/(2*alpha_d+5))*(r0**(2*alpha_d+5))

        ax6.scatter(ge,ann5)
        ax6.plot(ge,ge)
        ax6.plot(ge,0.1*ge)
        ax6.plot(ge,10*ge)
        axbonk(ax6,xlabel='GE',ylabel='all_density',xscale='log',yscale='log')







        fig.savefig('plots_to_sort/masses_%s'%sim)
        plt.close(fig)
        fig,ax=plt.subplots(1,1)
        ax.hist(agd,histtype='step')
        fig.savefig('plots_to_sort/a_from_b_hist_%s'%sim)
