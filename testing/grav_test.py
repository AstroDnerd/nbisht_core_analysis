from starter2 import *
import xtra_energy
import pcolormesh_helper as pch
import three_loopers_u500 as TL
sim_list = ['u501']
import gravity
for sim in sim_list:
    frame=10
    ds = TL.loops[sim].load(frame)
    xtra_energy.add_grav_test(ds)

    ad = ds.all_data()
    rho = ad[YT_density]
    arho = ad[YT_also_rho]
    #arho += rho.min()-arho.min()
    ext = extents()
    ext(rho)
    #ext(arho)
    #print(ext)

    if 0:
        proj = yt.ProjectionPlot(ds,0,YT_also_rho)
        proj.save('plots_to_sort/t1_%s_n%04d'%(sim,frame))
        proj = yt.ProjectionPlot(ds,0,YT_density)
        proj.save('plots_to_sort/t2_%s_n%04d'%(sim,frame))

    if 0:
        bins = np.geomspace( ext.minmax[0], ext.minmax[1],64)
        hist,xbins,ybins=np.histogram2d( rho, arho, bins=[bins,bins])
        fig,ax=plt.subplots(1,1)
        pch.helper(hist,xbins,ybins,ax=ax)
        ax.set(xscale='log',yscale='log',xlabel='rho',ylabel='grad g')
        ax.plot( ext.minmax,ext.minmax)
        fig.savefig('plots_to_sort/phase_%s_n%04dpng'%(sim,frame))

    if 0:
        cg = ds.covering_grid(0,[0.0]*3,[128]*3)
        rho = cg[YT_density]
        phi = cg[YT_potential_field]
        grav = gravity.gravity(rho.v, ds.parameters['GravitationalConstant'])
        #grav = gravity.gravity(rho.v, 1)
        grav.solve()
        ext=extents()
        phi1=phi.v.flatten()
        phi2=grav.phi.flatten()
        ext(phi1)
        ext(phi2)
        bins=np.linspace(ext.minmax[0],ext.minmax[1],64)

        hist,xbins,ybins=np.histogram2d( phi1, phi2, bins=[bins,bins])
        fig,ax=plt.subplots(1,1)
        pch.helper(hist,xbins,ybins,ax=ax)
        #ax.set(xscale='log',yscale='log',xlabel='phi disk',ylabel='phi fft')
        ax.plot( ext.minmax,ext.minmax)
        fig.savefig('plots_to_sort/phi_vs_phi_%s_n%04dpng'%(sim,frame))

    if 0:
        cg = ds.covering_grid(0,[0.0]*3,[128]*3)
        #rho = cg[YT_density]
        phi = cg[YT_potential_field]
        gdisk = cg[YT_acceleration_x]
        dx=1/128
        dphidx = -((phi[2:,:,:]-phi[:-2,:,:])/2)/dx
        print(dphidx.shape,"dphidx")
        g1 = dphidx[:,1:-1,1:-1]
        g2 = gdisk[1:-1,1:-1,1:-1]
        print(g1.shape)
        print(g2.shape)
        #a1=g1[0:2,0:2,0:2]
        #a2=g2[0:2,0:2,0:2]
        #print(a2/a1/(-colors.G))

        phi1=g1.flatten()
        phi2=g2.flatten()
        ext=extents()
        ext(phi1)
        reload(pch)
        fig,ax=plt.subplots(1,1)
        pch.simple_phase(phi1,phi2,ax=ax)
        ax.plot( ext.minmax,ext.minmax)
        fig.savefig('plots_to_sort/gx_vs_gx_b_%s_n%04dpng'%(sim,frame))

    if 0:
        #old plot method
        ext=extents()
        ext(phi1)
        ext(phi2)
        bins=np.linspace(ext.minmax[0],ext.minmax[1],64)

        hist,xbins,ybins=np.histogram2d( phi1, phi2, bins=[bins,bins])
        fig,ax=plt.subplots(1,1)
        pch.helper(hist,xbins,ybins,ax=ax)
        #ax.set(xscale='log',yscale='log',xlabel='phi disk',ylabel='phi fft')
        ax.plot( ext.minmax,ext.minmax)
        fig.savefig('plots_to_sort/gx_vs_gx_%s_n%04dpng'%(sim,frame))


    if 1:
        cg = ds.covering_grid(0,[0.0]*3,[128]*3)
        rho = cg[YT_density]
        phi = cg[YT_potential_field]
        gx = cg[YT_acceleration_x]
        gy = cg[YT_acceleration_y]
        gz = cg[YT_acceleration_z]
        dx=1/128
        dgxdx= ( gx[2:,:,:]-gx[:-2,:,:])/(2*dx)
        dgydy= ( gy[:,2:,:]-gy[:,:-2,:])/(2*dx)
        dgzdz= ( gz[:,:,2:]-gz[:,:,:-2])/(2*dx)

        f1 = 4*np.pi*colors.G*(rho[1:-1,1:-1,1:-1])
        f2 = -(dgxdx[:,1:-1,1:-1]+dgydy[1:-1,:,1:-1]+dgzdz[1:-1,1:-1,:])

    if 1:
        fig4,ax4=plt.subplots(1,2)
        ax4[0].imshow(f1.sum(axis=0))
        ax4[1].imshow(f2.sum(axis=0))
        fig4.savefig('plots_to_sort/proj')
    if 1:
        f1=f1.flatten()
        f2=f2.flatten()
        ext=extents()
        ext(f1)
        ext(f2)
        fig,ax=plt.subplots(1,1)
        pch.simple_phase(f1,f2,ax=ax)
        ax.plot( ext.minmax,ext.minmax)

        fig.savefig('plots_to_sort/n_vs_gx_%s_n%04dpng'%(sim,frame))

