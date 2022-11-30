from starter2 import *

import three_loopers_u500 as TL
import mountain_top

class br_cumsum():
    def __init__(self,this_looper):
        self.this_looper=this_looper
    def run(self, core_list=None,frame_list=None):
        this_looper=self.this_looper
        thtr=this_looper.tr


        if core_list is None:
            core_list=np.unique(this_looper.tr.core_ids)
        if frame_list is None:
            frame_list=this_looper.tr.frames

        for iframe, frame in enumerate(frame_list):
            ds = this_looper.load(frame)
            radius=1e-2
            radius = ds.arr(radius,'code_length')
            nframe=np.where(thtr.frames==frame)[0][0]

            for core_id in core_list:
                print('yay',core_id)
                ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True)

                for los in [1]:
                    normal = np.zeros(3)
                    normal[los]=1
                    XYZ = 'xyz'[los]
                    BI = 'magnetic_field_'+XYZ
                    p = np.array([ms.mean_x[nframe],ms.mean_y[nframe],ms.mean_z[nframe]])
                    c = ds.arr(p,'code_length')
                    R1 = 1./128
                    radius = ds.arr(R1,'code_length')
                    disk_long = ds.disk(c, normal, radius, 1.0)
                    disk_short = ds.disk(c, normal, radius, radius)

                    fig,axes=plt.subplots(2,2)
                    ax0=axes[0][0]; ax1=axes[1][0]; ax2=axes[0][1];ax3=axes[1][1]
                    arg_short = np.argsort( disk_short['enzo',XYZ])
                    arg_long = np.argsort( disk_long['enzo',XYZ])
                    rho_short = disk_short['density'][arg_short]
                    rho_long = disk_long['density'][arg_long]
                    z_short = disk_short['enzo',XYZ][arg_short]
                    z_long = disk_long['enzo',XYZ][arg_long]
                    B_short = disk_short[BI][arg_short]*rho_short
                    B_long = disk_long[BI][arg_long]*rho_long

                    mask = ms.compute_unique_mask(core_id, 1/2048, nframe)
                    Bi_part_only = thtr.c([core_id], BI)[mask,nframe]
                    rho_part= thtr.c([core_id], 'density')[mask,nframe]
                    dv      = thtr.c([core_id], 'cell_volume')[mask,nframe]
                    xyz_part = thtr.c([core_id], XYZ)[mask,nframe]
                    psort = np.argsort( xyz_part)
                    Bi_part = Bi_part_only*rho_part
                    rho_part = rho_part[psort]
                    Bi_part = Bi_part[psort]





                    crl=np.cumsum(rho_long)
                    crs=np.cumsum(rho_short)
                    cbl=np.cumsum( B_long)
                    cbs=np.cumsum( B_short)

                    start_index = np.argmin( np.abs( z_short.v-xyz_part.min()))
                    rho_off=crs[start_index].v
                    B_off = cbs[start_index].v


                    crp = np.cumsum( rho_part) + rho_off
                    cbp = np.cumsum( Bi_part) + B_off


                    ax0.plot( z_long,  crl, 'k:')
                    ax0.plot( z_short, crs,'k')
                    ax2.plot( z_long,  cbl,'k:')
                    ax2.plot( z_short, cbs,'k')
                    ax1.plot( z_long,  crl, 'k:')
                    ax1.plot( z_short, crs,'k')
                    ax3.plot( z_long,  cbl,'k:')
                    ax3.plot( z_short, cbs,'k')

                    ax0.plot( xyz_part[psort], crp, 'r')
                    ax2.plot( xyz_part[psort], cbp, 'r')



                    ax0.set(xlim=extents(z_short).minmax)
                    ax2.set(xlim=extents(z_short).minmax)
                    ax0.set(xlabel=XYZ,ylabel=r'$\int_0^%s \rho d%s$'%(XYZ,XYZ))
                    ax1.set(xlabel=XYZ,ylabel=r'$\int_0^%s \rho d%s$'%(XYZ,XYZ))
                    ax2.set(xlabel=XYZ,ylabel=r'$\int_0^%s \rho B_%s d%s$'%(XYZ,XYZ,XYZ))
                    ax3.set(xlabel=XYZ,ylabel=r'$\int_0^%s \rho B_%s d%s$'%(XYZ,XYZ,XYZ))
                    ax0.set_title(r'$\rho_s/\rho_l = %0.2f$'%(crs[-1]/crl[-1]))
                    ax2.set_title(r'$B_s/B_l = %0.2f$'%(cbs[-1]/cbl[-1]))
                    fig.tight_layout()
                    fig.savefig('plots_to_sort/rays_c%04d_%s'%(core_id, XYZ))

class br_mean_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
        self.abs_b_long=[]
        self.abs_b_short=[]
        self.mean_b_long=[]
        self.mean_b_short=[]
        self.mean_b_part=[]
        self.mean_rho_long=[]
        self.mean_rho_short=[]
        self.mean_rho_part=[]
        self.long_volume=[]
        self.short_volume=[]
    def run(self, core_list=None,frame_list=None):
        this_looper=self.this_looper
        thtr=this_looper.tr


        if core_list is None:
            core_list=np.unique(this_looper.tr.core_ids)
        if frame_list is None:
            frame_list=this_looper.tr.frames

        for iframe, frame in enumerate(frame_list):
            ds = this_looper.load(frame)
            radius=1e-2
            radius = ds.arr(radius,'code_length')
            nframe=np.where(thtr.frames==frame)[0][0]

            for core_id in core_list:
                print('yay',core_id)
                ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True)

                for los in [1]:
                    normal = np.zeros(3)
                    normal[los]=1

                    XYZ = 'xyz'[los]
                    BI = 'B'+XYZ



                    p = np.array([ms.mean_x[nframe],ms.mean_y[nframe],ms.mean_z[nframe]])
                    c = ds.arr(p,'code_length')
                    R1 = 1./128
                    radius = ds.arr(R1,'code_length')
                    disk_long = ds.disk(c, normal, radius, 1.0)
                    disk_short = ds.disk(c, normal, radius, radius)

                    mass_or_volume='gas','cell_volume'
                    total_long = disk_long[mass_or_volume].sum()
                    total_short = disk_short[mass_or_volume].sum()
                    self.long_volume.append(total_long)
                    self.short_volume.append(total_short)
                    total_long=1.0
                    total_short=1.0
                    self.mean_b_long.append( (disk_long[BI]*disk_long['cell_mass']).sum()/disk_long['cell_mass'].sum())
                    self.mean_b_short.append( (disk_short[BI]*disk_short['cell_mass']).sum()/disk_short['cell_mass'].sum())
                    self.mean_rho_long.append(   (disk_long['density']*disk_long[mass_or_volume]).sum()/total_long) 
                    self.mean_rho_short.append( (disk_short['density']*disk_short[mass_or_volume]).sum()/total_short) 

                    #particles
                    Bi_part = thtr.c([core_id], BI)[:,nframe]
                    rho_part= thtr.c([core_id], 'density')[:,nframe]
                    dv      = thtr.c([core_id], 'cell_volume')[:,nframe]
                    total_vol = dv.sum()
                    total_mass = (rho_part*dv).sum()
                    self.mean_rho_part.append(  (rho_part*dv).sum()/total_vol)
                    self.mean_b_part.append( (Bi_part*rho_part*dv).sum()/total_mass)


                    if 1:
                        fig,axes=plt.subplots(2,2)
                        ax0=axes[0][0]; ax1=axes[1][0]; ax2=axes[0][1];ax3=axes[1][1]
                        arg_short = np.argsort( disk_short['enzo',XYZ])
                        arg_long = np.argsort( disk_long['enzo',XYZ])
                        rho_short = disk_short['density'][arg_short]
                        rho_long = disk_long['density'][arg_long]
                        z_short = disk_short['enzo',XYZ][arg_short]
                        z_long = disk_long['enzo',XYZ][arg_long]
                        B_short = disk_short['Bx'][arg_short]*rho_short
                        B_long = disk_long['Bx'][arg_long]*rho_long
                        ax0.plot( z_long, np.cumsum(rho_long), 'k:')
                        ax0.plot( z_short, np.cumsum(rho_short),'k')
                        ax2.plot( z_long, np.cumsum( B_long)/np.cumsum(rho_long),'k:')
                        ax2.plot( z_short, np.cumsum( B_short)/np.cumsum(rho_short),'k')
                        ax0.set(xlim=extents(z_short).minmax)
                        ax2.set(xlim=extents(z_short).minmax)
                        fig.savefig('plots_to_sort/rays_c%04d'%core_id)
                    #proj = ds.proj('density',los,center=c,data_source=disk)
                    #dx=1./2048
                    #nx = np.ceil(2*radius/dx)
                    #frb = proj.to_frb(2*radius,nx, center=c)
                    #fig,axS=plt.subplots(2,2)
                    #ax =axS.flatten()
                    #plot=ax[0].imshow(np.log10(frb['density']))
                    #fig.colorbar(plot,ax=ax[0])
                    #ax[0].set_title('rho')
                    #norm=mpl.colors.LogNorm(vmin=bmin,vmax=bmax)
                    #Btotal=frb['magnetic_field_strength']
                    #b_ext(Btotal[Btotal>0])
                    #plot=ax[1].imshow(Btotal)#,norm=norm)
                    #fig.colorbar(plot,ax=ax[1])
                    #ax[1].set_title('Btotal')
                    #ax[2].imshow((frb['magnetic_field_z']))

                    #xbins = np.linspace( pos[0].min(),pos[0].max(),129)
                    #ybins = np.linspace( pos[1].min(),pos[1].max(),129)
                    #hist, xb, yb = np.histogram2d( pos[0].flatten(), pos[1].flatten(), bins=[xbins,ybins],weights=dv.flatten())
                    ##pch.helper(hist.transpose(), yb.transpose(), xb.transpose(), ax=ax0)
                    #pch.helper(hist, xb, yb, ax=ax0, transpose=False)
            
           


                    #fig.savefig('plots_to_sort/B-rho_%s_c%04d_n%04d.png'%(this_looper.sim_name,core_id,frame))
                    #pw = proj.to_pw()
                    #pw.zoom(0.5/radius.v)
                    #pw.save('plots_to_sort/thing_%s_c%04d'%(this_looper.sim_name,core_id))


if 1:
    #cumsums
    for nsim,sim in enumerate(TL.loops):
        if nsim != 0:
            continue
        this_looper = TL.loops[sim]
        core_list=np.unique(this_looper.tr.core_ids)
        core_list=core_list[15:20]
        #core_list=[323]
        #core_list=None
        frames=this_looper.tr.frames[-1:]
        #proj=brho(this_looper,core_list=core_list,frame_list=frames)
        BRC = br_cumsum(this_looper)
        BRC.run(core_list=core_list, frame_list=frames)

if 0:
    if 'BRT' not in dir():
        for nsim,sim in enumerate(TL.loops):
            if nsim != 0:
                continue
            this_looper = TL.loops[sim]
            core_list=np.unique(this_looper.tr.core_ids)
            #core_list=core_list[:20]
            #core_list=[323]
            #core_list=None
            frames=this_looper.tr.frames[-1:]
            #proj=brho(this_looper,core_list=core_list,frame_list=frames)
            BRT = br_mean_tool(this_looper)
            BRT.run(core_list=core_list, frame_list=frames)

    if 1:
        fig,axes=plt.subplots(2,2)
        ax0=axes[0][0];ax1=axes[0][1];ax2=axes[1][0];ax3=axes[1][1]
        rho_L = nar(BRT.mean_rho_long)
        rho_S = nar(BRT.mean_rho_short)
        B_L = nar(BRT.mean_b_long)
        B_S = nar(BRT.mean_b_short)
        VL = nar(BRT.long_volume)
        VS = nar(BRT.short_volume)

        ax0.scatter( rho_S,rho_L)
        xtr = extents(rho_S)
        xtr(rho_L)
        ax0.plot(xtr.minmax,xtr.minmax)
        ax0.set(xlabel='rhoS',ylabel='rhoL')

        ax1.scatter( B_S,B_L)
        ax1.set(xlabel='BS',ylabel='BL')
        xtB = extents(B_S)
        xtB(B_L)
        ax1.plot(xtB.minmax,xtB.minmax)

        ax2.scatter( rho_S/rho_L,B_S/B_L )
        ax2.set(xlabel='rhoS/rhoL',ylabel='BS/BL')
        ax2.axhline(1, c=[0.5]*4)


        Q1 = B_S/B_L
        Q2 = rho_S/rho_L
        V = np.linspace(0,1,Q1.size)
        ax3.plot( sorted(Q1), V)
        ax3.plot( sorted(Q2), V)
        ax3.axvline(1,c=[0.5]*3)
        ax3.set(xlabel='Q_S/Q_L',ylabel='N')
        #ax3.legend(loc=0)

        #xtB = extents(B_S)
        #xtB(B_L)
        #ax3.plot(xtB.minmax,xtB.minmax)

        #ax2.scatter( VS,VL)
        #ax2.set(xlabel='VS',ylabel='VL')
        #xt = extents()
        #xt(X)
        #xt(Y)
        #ax2.plot( xt.minmax,xt.minmax)
        #X,Y=rho_S/rho_L, B_S/B_L
        #ax2.scatter( X,Y)
        fig.tight_layout()
        fig.savefig('plots_to_sort/test.png')

