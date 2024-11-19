from starter2 import *
import track_loader as TL
#import ucore
#reload(ucore)

import these_particles5


import dtools.davetools as dt
import pcolormesh_helper as pch


sim_list=['u501']#,'u502','u503']
import track_loader as TL
TL.load_tracks(sim_list)
import monster
monster.load(sim_list)

def add_sphere(axes,direction,radius,center,color):
    xx = [1,0,0][direction]
    yy = [2,2,1][direction]
    cen = nar([center[xx],center[yy]])
    circle = plt.Circle(cen,radius,fill=False,color=color)
    axes.add_artist(circle)

def add_particles(axes,direction,positions,color):
    xx = [1,0,0][direction]
    yy = [2,2,1][direction]
    px = positions[xx]
    py = positions[yy]
    axes.scatter(px,py,c=color,s=1,alpha=0.5)




def image(mon,core_list,frames):

    if type(frames) == str:
        if frames == 'movie':
            frame_slice=slice(None)
        elif frames == 'short':
            frame_slice = np.zeros_like(mon.frames,dtype='bool')
            frame_slice[::10]=True
            frame_slice[-1]=True
        elif frames == 'realshort':
            frame_slice=slice(-1,None,10)
            frame_slice = np.zeros_like(mon.frames,dtype='bool')
            frame_slice[0]=True;frame_slice[-1]=True
            #frame_slice[-2]=True
        frame_list=mon.frames[frame_slice]


    for frame in frame_list:
        for core_id in core_list:
            print('image',mon.name,core_id,frame)
            sph1 = mon.get_sphere(core_id,frame,'r1')
            sphmax = mon.get_sphere(core_id,frame,'rmax')
            sphinf = mon.get_sphere(core_id,frame,'rinf')
            sphsmrt = mon.get_sphere(core_id,frame,'rsmart')
            rad = [sph1.radius, sphmax.radius,sphinf.radius]
            sss = [sph1, sphmax,sphinf]
            which = np.argmax(rad)
            sph = sss[which]
            ds = mon.get_ds(frame)
            direction=2
            proj = ds.proj('density',direction,data_source=sph, center=sph.center)
            frb = proj.to_frb(2*sph.radius,512)
            rho=frb['density']
            left = sph.center-sph.radius
            right = sph.center+sph.radius
            xx = [1,0,0][direction]
            yy = [2,2,1][direction]

            cmap=copy.copy(mpl.cm.get_cmap("gray"))
            cmap.set_under('w')
            norm = mpl.colors.LogNorm(vmin=rho[rho>0].min(),vmax=rho.max())

            fig,ax=plt.subplots(1,2,figsize=(8,4))
            ax0=ax[0];ax1=ax[1]#;ax2=ax[2]

            ax0.imshow(frb['density'].transpose(),extent=[left[xx],right[xx],left[yy],right[yy]],norm=norm,cmap=cmap)
            add_sphere(ax0,direction,sph1.radius,sph.center,'r')
            add_sphere(ax0,direction,sphmax.radius,sph.center,'g')
            add_sphere(ax0,direction,sphinf.radius,sph.center,'orange')

            ms = mon.get_ms(core_id)
            nf = mon.get_frame_index(frame)
            positions = np.stack([ms.this_x[:,nf],ms.this_y[:,nf],ms.this_z[:,nf]])
            add_particles(ax0,direction,positions,'purple')

            center = sph.center+0
            center.shape = center.size,1
            re_cen = positions-center.v
            r = ((re_cen**2).sum(axis=0))**0.5
            r[r<1/2048]=1/2048
            
            #the_r_p = ms.r[:,nf]
            the_r_p = r
            the_d_p = ms.density[:,nf]
            the_r_g = sph['radius']
            the_d_g = sph['density']
            ext_r = dt.extents()
            ext_r(the_r_p)
            ext_r(the_r_g)
            ext_d = dt.extents()
            ext_d(the_d_p)
            ext_d(the_d_g)
            rbins = np.geomspace(1/2048, ext_r.minmax[1],64)
            dbins = np.geomspace(ext_d.minmax[0],ext_d.minmax[1],64)
            pch.simple_phase( sph['radius'],sph['density'],bins=[rbins,dbins],ax=ax1)
            ax1.set(xscale='log',yscale='log')
            
            ax1.scatter(the_r_p,the_d_p,c='r')
            ax1.axvline(sph1.radius,c='r')
            ax1.axvline(sphmax.radius,c='g')
            ax1.axvline(sphinf.radius,c='orange')

            ax1b=ax1.twinx()
            the_x = r+0
            ORDER = np.argsort(the_x)
            dv = ms.cell_volume[:,nf][ORDER]
            rho_srt = the_d_p[ORDER]
            m = np.cumsum(dv*rho_srt)
            the_y = np.arange(1,r.size+1)/r.size
            #ax1b.plot(the_x,the_y)
            ax1b.plot(the_x[ORDER],m)

            fig.savefig("%s/proj_%s_c%04d_n%04d"%(plot_dir,mon.name,core_id,frame))




for sim in sim_list[:1]:
    mon = monster.closet[sim]
    core_list =  mon.this_looper.core_by_mode['A'][3:4]
    pw=image(mon,core_list,frames='short')
