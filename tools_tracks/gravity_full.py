from starter2 import *

import three_loopers_six as TL
simlist=['u603']
core_id = 1
import gravity
reload(gravity)

import pcolormesh_helper as pch
reload(pch)
plt.close('all')

extract_level=2; extract_size=512
#extract_level=3; extract_size=1024
t0=time.time()
print('start')
for sim in simlist:
    frame = TL.loops[sim].target_frame
    this_looper = TL.loops[sim]
    ds = TL.loops[sim].load(frame)
    print('make cg')
    if 'cg' not in dir():
        cg = ds.covering_grid(extract_level,[0.0]*3,[extract_size]*3)
        rho=cg[YT_density]
    else:
        print('use old cg')
    if 'ggg' not in dir():
        ggg = gravity.gravity( rho, ds.parameters['GravitationalConstant'])
        print('do solve')
        ggg.solve()
        print('done')
    else:
        print("WARNING old ggg")

    
    if 1:
        ms = trackage.mini_scrubber( this_looper.tr, core_id)
        c = nar([ms.mean_x[-1], ms.mean_y[-1],ms.mean_z[-1]])


        if 'rho_hack' not in dir():
            print("On board the memory train") 
            rho_hack = rho+0
            x = cg[YT_x].v; y=cg[YT_y].v; z=cg[YT_z].v
            r = np.sqrt((x-c[0])**2 + (y-c[1])**2 + (z-c[2])**2)
            density_cut = rho > 1e4
            radius_cut = r<0.01
            rho_hack[ density_cut*radius_cut]=0
            print("Success!")
        else:
            print('old rho_hack')

        if 'ggg2' not in dir():
            print('do new solve')
            ggg2=gravity.gravity( rho_hack, ds.parameters['GravitationalConstant'])
            ggg2.solve()
            print('done')
        else:
            print("OLD GGG2")

    if 1:
        print('plotting')
        ii = (c*extract_size).astype('int')

        fig,ax=plt.subplots(3,2,figsize=(12,12))
        ax0=ax[0][0]; ax1=ax[0][1]; ax2=ax[1][0]
        ax3=ax[1][1]; ax4=ax[2][0]; ax5=ax[2][1]

        def fix(q):
            qout = max([q,0])
            qout = min([qout, extract_size])
            return qout
        ix=ii[0]
        ix0=fix(ii[0]-20)
        ix1=fix(ii[0]+20)
        iy0=fix(ii[1]-20)
        iy1=fix(ii[1]+20)
        iz0=fix(ii[2]-20)
        iz1=fix(ii[2]+20)

        print('proj rho')
        S1 = cg[YT_density][ix,iy0:iy1,iz0:iz1]
        norm = mpl.colors.LogNorm( vmin=S1.min(), vmax=S1.max())
        ploot=ax0.imshow(S1,norm=norm)
        fig.colorbar(ploot,ax=ax0)
        S2 = rho_hack[ix,iy0:iy1,iz0:iz1]
        ploot=ax1.imshow(S2,norm=norm)
        fig.colorbar(ploot,ax=ax1)


        if 0:
            S3 = ggg.phi[ix,iy0:iy1,iz0:iz1]
            norm = mpl.colors.Normalize( vmin=S3.min(), vmax=S3.max())
            ploot=ax3.imshow(S3,norm=norm)
            fig.colorbar(ploot,ax=ax3)
            S2 = ggg2.phi[ix,iy0:iy1,iz0:iz1]
            ploot=ax4.imshow(S2,norm=norm)
            fig.colorbar(ploot,ax=ax4)


        print('make things')
        phi1 = ggg.phi[ix0:ix1,iy0:iy1,iz0:iz1]
        phi2 = ggg2.phi[ix0:ix1,iy0:iy1,iz0:iz1]
        phi_min,phi_max= min([phi1.min(),phi2.min()]),max([phi1.max(),phi2.max()])
        rcube=r[ix0:ix1,iy0:iy1,iz0:iz1]
        dv_cube = cg[YT_cell_volume][ix0:ix1,iy0:iy1,iz0:iz1]

        ps1= phi1.sum(axis=0)
        ps2= phi2.sum(axis=0)
        p3min,p3max= min([ps1.min(),ps2.min()]),max([ps1.max(),ps2.max()])
        norm = mpl.colors.Normalize(p3min,p3max)
        ax2.imshow(ps1,norm=norm)
        ax3.imshow(ps2,norm=norm)

        print('more plots')
        if 1:
            rbins = np.geomspace(1/2048, rcube.max(), 64)
            phibins=np.linspace( phi_min, phi_max, 65)
            print('  hist1')
            hist1, xb1, yb1 = np.histogram2d( rcube.flatten(), phi1.flatten(),
                                             bins=[rbins,phibins], weights=dv_cube.flatten())
            pch.helper(hist1,xb1,yb1,ax=ax4)
            print('  hist2')
            hist2, xb2, yb2 = np.histogram2d( rcube.flatten(), phi2.flatten(), 
                                             bins=[rbins,phibins], weights=dv_cube.flatten())
            pch.helper(hist2,xb2,yb2,ax=ax5)
            ax4.set_xscale('log')
            ax5.set_xscale('log')
            ax4.set_ylim(phi_min,phi_max)
            ax5.set_ylim(phi_min,phi_max)


        if 0:
            ax4.scatter( rcube.flatten(), phi1.flatten(),c='k',s=0.1)
            ax5.scatter( rcube.flatten(), phi2.flatten(),c='r',s=0.1)
            #ax4.set_ylim(phi_min,phi_max)
            #ax5.set_ylim(phi_min,phi_max)
            ax4.set_xscale('log')
            ax5.set_xscale('log')


        if 0:
            ppmin=min([phi1.min(),phi2.min()])
            ppmax=max([phi1.max(),phi2.max()])
            phibins = np.linspace( ppmin,ppmax,64)
            hist, xb, yb = np.histogram2d( phi1, phi2, bins=[phibins,phibins], weights=dv)
            pch.helper(hist,xb,yb,ax=ax5)


        print('savefig')
        fig.savefig('plots_to_sort/rho_%s_c%04d.png'%(sim,core_id))
        plt.close('fig')

t1=time.time()
print('TIME',t1-t0)
if 0:
    if 0:
        #projections of potential.
        fig,ax=plt.subplots(2,2,figsize=(12,4))
        ax0=ax[0][0]; ax1=ax[0][1]; ax2=ax[1][0]; ax3=ax[1][1]
        ploot=ax0.imshow( cg[YT_potential_field].sum(axis=0))
        fig.colorbar(ploot,ax=ax0)
        ploot=ax1.imshow( ggg.phi.sum(axis=0))
        fig.colorbar(ploot,ax=ax1)
        fig.savefig('plots_to_sort/grav_grav_%s.png'%this_looper.sim_name)
        plt.close(fig)


    if 0:
        fig,ax=plt.subplots(1,1)
        phi1 = cg[YT_potential_field].flatten().v
        phi2 = ggg.phi.flatten()
        #phi vs phi. Works pretty well.
        phibins = np.linspace( min([phi1.min(),phi2.min()]), max([phi1.max(),phi2.max()]),64)
        dv = cg[YT_cell_volume].v.flatten()
        hist, xb, yb = np.histogram2d( phi1, phi2, bins=[phibins,phibins], weights=dv)
        pch.helper(hist,xb,yb,ax=ax)
        #ax.scatter( cg[YT_density].flatten(), ggg.phi.flatten())
        fig.savefig('plots_to_sort/ggg.png')

    if 0:
        am = np.argmax(rho)
        x = cg[YT_x]; y=cg[YT_y]; z=cg[YT_z]

        ix = int(x.flatten()[am]*extract_size)
        fig,ax=plt.subplots(1,2,figsize=(12,4))

        S1 = cg[YT_potential_field][ix,:,:]
        ploot=ax[0].imshow(S1)
        fig.colorbar(ploot,ax=ax[0])
        S2 = ggg.phi[ix,:,:]
        ploot=ax[1].imshow(S2)
        fig.colorbar(ploot,ax=ax[1])
        fig.savefig('plots_to_sort/phi_max_%s.png'%sim)
        plt.close('fig')

    if 0:
        ms = trackage.mini_scrubber( this_looper.tr, core_id)
        c = nar([ms.mean_x[-1], ms.mean_y[-1],ms.mean_z[-1]])

        ii = (c*extract_size).astype('int')

        fig,ax=plt.subplots(1,3,figsize=(12,4))

        ix=ii[0]
        iy0=ii[1]-20
        iy1=ii[1]+20
        iz0=ii[2]-20
        iz1=ii[2]+20


        S1 = cg[YT_potential_field][ix,iy0:iy1,iz0:iz1]
        ploot=ax[0].imshow(S1)
        fig.colorbar(ploot,ax=ax[0])
        S2 = ggg.phi[ix,iy0:iy1,iz0:iz1]
        ploot=ax[1].imshow(S2)
        fig.colorbar(ploot,ax=ax[1])


        ax[2].scatter( S1, S2)
        ax[2].plot(S1.flatten(),S1.flatten())
        fig.savefig('plots_to_sort/rho_%s_c%04d.png'%(sim,core_id))
        plt.close('fig')
