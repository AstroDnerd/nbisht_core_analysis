from starter2 import *

import three_loopers_u500 as TL
import mountain_top

def brho(this_looper, core_list=None,frame_list=None):

    if core_list is None:
        core_list=np.unique(this_looper.tr.core_ids)
    if frame_list is None:
        frame_list=this_looper.tr.frames
    thtr=this_looper.tr

    reload(mountain_top)
    b_ext=extents()
    for iframe, frame in enumerate(frame_list):
        ds = this_looper.load(frame)
        radius=1e-2
        radius = ds.arr(radius,'code_length')
        nframe=np.where(thtr.frames==frame)[0][0]
        bx = thtr.track_dict['magnetic_field_x']
        by = thtr.track_dict['magnetic_field_y']
        bz = thtr.track_dict['magnetic_field_z']
        B = np.sqrt(bx**2+by**2+bz**2)
        bmin=B[B>0].min()
        bmax=B.max()

        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            bx = thtr.c([core_id],'magnetic_field_x')
            by = thtr.c([core_id],'magnetic_field_y')
            bz = thtr.c([core_id],'magnetic_field_z')
            B = np.sqrt(bx**2+by**2+bz**2)
            bmin=B[B>0].min()
            bmax=(B.max())**0.75
            bmin=10
            bmax=20

            for los in [0]:
                normal = np.zeros(3)
                normal[los]=1
                p = np.array([ms.mean_x[nframe],ms.mean_y[nframe],ms.mean_z[nframe]])
                c = ds.arr(p,'code_length')
                radius = ds.arr(5e-2,'code_length')
                if 0:
                    disk = ds.disk(c, normal, radius, 1)
                elif 0:
                    L = c-radius
                    R = c+radius
                    L[los]=0
                    R[los]=1
                    disk=ds.region(c,L,R)
                else:
                    disk = ds.sphere(c,radius)

                proj = ds.proj('density',los,center=c,data_source=disk)
                dx=1./2048
                nx = np.ceil(2*radius/dx)
                frb = proj.to_frb(2*radius,nx, center=c)
                fig,axS=plt.subplots(2,2)
                ax =axS.flatten()
                plot=ax[0].imshow(np.log10(frb['density']))
                fig.colorbar(plot,ax=ax[0])
                ax[0].set_title('rho')
                norm=mpl.colors.LogNorm(vmin=bmin,vmax=bmax)
                Btotal=frb['magnetic_field_strength']
                b_ext(Btotal[Btotal>0])
                plot=ax[1].imshow(Btotal)#,norm=norm)
                fig.colorbar(plot,ax=ax[1])
                ax[1].set_title('Btotal')
                ax[2].imshow((frb['magnetic_field_z']))

                #xbins = np.linspace( pos[0].min(),pos[0].max(),129)
                #ybins = np.linspace( pos[1].min(),pos[1].max(),129)
                #hist, xb, yb = np.histogram2d( pos[0].flatten(), pos[1].flatten(), bins=[xbins,ybins],weights=dv.flatten())
                ##pch.helper(hist.transpose(), yb.transpose(), xb.transpose(), ax=ax0)
                #pch.helper(hist, xb, yb, ax=ax0, transpose=False)
        
       


                fig.savefig('plots_to_sort/B-rho_%s_c%04d_n%04d.png'%(this_looper.sim_name,core_id,frame))
                #pw = proj.to_pw()
                #pw.zoom(0.5/radius.v)
                #pw.save('plots_to_sort/thing_%s_c%04d'%(this_looper.sim_name,core_id))
    print(b_ext)

if 1:
    for nsim,sim in enumerate(TL.loops):
        if nsim != 0:
            continue
        this_looper = TL.loops[sim]
        core_list=np.unique(this_looper.tr.core_ids)
        core_list=core_list[5:6]
        core_list=[323]
        frames=this_looper.tr.frames[-30:]
        proj=brho(this_looper,core_list=core_list,frame_list=frames)



    #this_looper.targets

    #plot_peaks(this_looper, core_list=[323])

