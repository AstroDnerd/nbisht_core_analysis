from starter2 import *
import xtra_energy
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
sim_list=['u501']

from collections import defaultdict
import r_inflection

import three_loopers_u500 as TL
#import three_loopers_six as TL
if 'inflection' not in dir():
    inflection = {}
    for sim in TL.loops:
        inflection[sim]=r_inflection.R_INFLECTION( TL.loops[sim])

class slope_tool():
    def __init__(self, this_looper, inflection):
        self.this_looper=this_looper
        self.output = defaultdict(list)
        self.rinflection=inflection.rinflection
        self.rinflection_list=inflection.rinflection_list

    def run(self,core_list=None, do_plots=True,do_proj=True):
        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        this_looper=self.this_looper
        frame = this_looper.target_frame
        ds = this_looper.load(frame)
        G = ds['GravitationalConstant']/(4*np.pi)
        xtra_energy.add_energies(ds)
        for core_id in core_list:
            print('Potential %s %d'%(this_looper.sim_name,core_id))

            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            c = nar([ms.mean_x[-1], ms.mean_y[-1],ms.mean_z[-1]])
            
            R_SPHERE = 8/128
            rsph = ds.arr(R_SPHERE,'code_length')
            sp = ds.sphere(c,rsph)
            print('V1')




            GE = np.abs(sp['grav_energy'])
            dv = np.abs(sp['cell_volume'])
            RR = sp['radius']
            DD = sp['density']

            R_KEEP = R_SPHERE #self.rinflection[core_id]


            #get the zones in the sphere that are
            #within R_KEEP
            ok_fit = np.logical_and(RR  < R_KEEP, GE>0)

            rok=RR[ok_fit].v
            print('V2')

            #
            # Binding energy and GMM/R
            #
            ge_total = (GE[ok_fit]*dv[ok_fit]).sum()
            mtotal = (sp['cell_mass'][ok_fit]).sum()
            gmm      = G*mtotal**2/R_KEEP
            self.output['ge_total'].append(ge_total)
            self.output['gmm'].append(gmm)

            #store stuff
            self.output['mass'].append(mtotal)
            self.output['r0'].append(R_KEEP)
            print('V3')

            fig2,ax2=plt.subplots(1,1)
            #ay0=ax2[0][0]; ay1=ax2[0][1];ay2=ax2[1][0]ay3=ax3[1][1]
            ay0=ax2


            ORDER = np.argsort(rok)
            rho_o = DD[ok_fit][ORDER].v
            ge_o = GE[ok_fit][ORDER].v
            dv_o  = dv[ok_fit][ORDER].v
            rr_o  = RR[ok_fit][ORDER].v
            mass_r = (rho_o*dv_o).cumsum()
            enrg_r = (ge_o*dv_o).cumsum()
            print('V4')

            ay0.plot( rr_o, enrg_r, marker='*', c='k')

            #ge_r= 4*np.pi*A/( (rr_o**(2*alpha+2)*(2*alpha+5)))*(rr_o**(2*alpha+5)-rmin**(2*alpha+5))
            #ay0.plot( rr_o, ge_r)

            #ay0.scatter( R_KEEP, ge_good, c='r', marker='*')

            gmm_r = G*mass_r**2/rr_o
            ay0.scatter(R_KEEP,self.output['gmm'][-1])
            ay0.plot( rr_o, gmm_r, c='r')

            
            az0=ay0.twinx()
            all_r,all_m=rr_o[1:], mass_r[1:]
            from scipy.interpolate import interp1d
            mfunc = interp1d( all_r, all_m)
            my_r = np.linspace(1/2048,all_r.max(),128)
            my_m = mfunc(my_r)
            #ay0.plot( my_r, my_m * enrg_r.max()/my_m.max(), c='m')
            az0.plot( my_r, my_m )
            #az0.plot( rr_o, rho_o)
            #az0.set_yscale('log')
            dm=(my_m[1:]- my_m[:-1])
            mm =0.5*(my_m[1:]+my_m[:-1])
            dr=(my_r[1:]-my_r[:-1])
            dm_dr = dm/dr
            rbins = 0.5*(my_r[1:]+my_r[:-1])

            fig3,ax3=plt.subplots(1,1)
            SWITCH=2*rbins*dm_dr/mm 
            ax3.plot( rbins, SWITCH)
            ax3.plot( rbins, rbins*0+1)
            ax4=ax3.twinx()
            ax4.plot( my_r, my_m)
            ok = np.where(np.logical_and( rbins>0.01, SWITCH < 1))[0][0]
            print('OK',ok)
            ax4.scatter( my_r[ok-1], my_m[ok-1])

            fig3.savefig('plots_to_sort/dr.png')


            print('V5')

            #fig2.savefig('plots_to_sort/%s_cuml_c%04d.png'%(this_looper.sim_name, core_id))




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

                xx = ds.coordinates.x_axis[proj_axis]
                yy = ds.coordinates.y_axis[proj_axis]
                vh=frb['velocity_%s'%'xyz'[xx]][::8]
                vv=frb['velocity_%s'%'xyz'[yy]][::8]
                nx,ny=column_density.shape
                the_x,the_y = np.mgrid[0:nx:8,0:ny:8]
                #pdb.set_trace()
                #the_x=frb['xyz'[xx]]
                #the_y=frb['xyz'[yy]]

                rx1.imshow(column_density,norm=norm)
                rx1.quiver(the_x,the_y,vh,vv)

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

                if 1:
                    #plot GE
                    #ax0.plot(r_cen[keepers], UE)
                    pch.helper(h2,xb,yb,ax=ax0,transpose=False)
                    #ax0.plot( rok, GE_fit_line,c='r')
                    ax0.plot( rr, GE_fit_line,c='r')
                    #ax0.scatter( R_KEEP,UE[index],c='r')
                    #AlsoA=colors.G*mtotal**2/(4*np.pi)*R_KEEP**(-4)
                    ax0.scatter( R_KEEP, A, c='r',marker='*')
                    #ax0.scatter( R_KEEP, AlsoA, c='b',marker='*')
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



                outname='plots_to_sort/%s_c%04d_potfit'%(this_looper.sim_name,core_id)
                axbonk(ax0,xscale='log',yscale='log',xlabel='r',ylabel='grad phi sq')
                fig.savefig(outname)
                print(outname)

if 'stuff' not in dir() or True:
    stuff={}
    for sim in ['u503']:
        #stuff={}
        all_cores=np.unique( TL.loops[sim].tr.core_ids)
        core_list=list(all_cores)
        core_list=all_cores[:1]
        #core_list=[323]
        #stuff[sim]=plot_phi( TL.loops[sim],core_list=core_list, do_plots=False)
        stuff[sim]=slope_tool(TL.loops[sim], inflection[sim])
        stuff[sim].run(core_list=core_list,do_plots=False, do_proj=False)

