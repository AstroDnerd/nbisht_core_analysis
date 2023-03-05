from starter2 import *

import three_loopers_u900 as TL

from scipy import ndimage

import movie_frames
reload(movie_frames)
import heat_map
reload(heat_map)
plt.close('all')

class dq_dt():
    def __init__(self,this_looper):
        self.this_looper=this_looper
    def run(self,core_list=None,frame_list=None):
        this_looper=self.this_looper
        thtr = this_looper.tr

        if frame_list is None:
            mask = movie_frames.quantized_mask(this_looper)
            times = thtr.times[mask]
            if times[0] == times[1]:
                mask[0]=False
            times = thtr.times[mask]
            frame_list=this_looper.tr.frames[mask]

            mask_m1 = mask[:-1]
            mask_p1 = mask[1:]
            tcenter = 0.5*(times[:-1]+times[1:])
            dt = times[1:]-times[:-1]

            dt_square = dt+0
            dt_square.shape = dt_square.size,1


        if core_list is None:
            core_list = sorted(np.unique( thtr.core_ids))


        for core_id in core_list:
            print('go',core_id)
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            ms.particle_pos(core_id)

            fig,ax_square=plt.subplots(2,2)
            ax=ax_square.flatten()

            if ms.nparticles < 1000:
                sl=slice(None)
                c=[0.5]*4
            else:
                sl = slice(None,None,10)
                #c=[0,0,0,0.1]
                c=[0.1]*4
            rho = ms.density[sl].transpose()
            rho = rho[mask,:]
            Bmag=thtr.c([core_id],'magnetic_field_strength')[sl].transpose()[mask,:]#/colors.mean_field[this_looper.sim_name]
            BP = B**2/2
            B2 = B**2
            divv=thtr.c([core_id],'velocity_divergence')[sl].transpose()[mask,:]#/colors.mean_field[this_looper.sim_name]

            fig, ax=plt.subplots(3,6,figsize=(20,12))
            ax0=ax[0][0];ax1=ax[0][1] 
            ax2=ax[1][0];ax3=ax[1][1]
            ax4=ax[2][0];ax5=ax[2][1]
            ax6=ax[0][2];
            ax7=ax[1][2]
            ax8=ax[2][2]
            ax9=ax[0][3]
            ax10=ax[1][3]
            ax11=ax[2][3]
            ax12=ax[0][4]
            ax13=ax[1][4]
            ax14=ax[2][4]
            ax15=ax[0][5]
            ax16=ax[1][5]
            ax17=ax[2][5]

            def extract(arr):
                return arr[sl].transpose()[mask,:]
            bx   =extract(thtr.c([core_id],'magnetic_field_x'))
            by   =extract(thtr.c([core_id],'magnetic_field_y'))
            bz   =extract(thtr.c([core_id],'magnetic_field_z'))
            vx   =extract(thtr.c([core_id],'velocity_x'))
            vy   =extract(thtr.c([core_id],'velocity_y'))
            vz   =extract(thtr.c([core_id],'velocity_z'))
            dxvx =extract(thtr.c([core_id],'dxvx'))
            dxvy =extract(thtr.c([core_id],'dxvy'))
            dxvz =extract(thtr.c([core_id],'dxvz'))
            dyvx =extract(thtr.c([core_id],'dyvx'))
            dyvy =extract(thtr.c([core_id],'dyvy'))
            dyvz =extract(thtr.c([core_id],'dyvz'))
            dzvx =extract(thtr.c([core_id],'dzvx'))
            dzvy =extract(thtr.c([core_id],'dzvy'))
            dzvz =extract(thtr.c([core_id],'dzvz'))
            Sx = bx*dxvx+by*dyvx+bz*dzvx
            Sy = bx*dxvy+by*dyvy+bz*dzvy
            Sz = bx*dxvz+by*dyvz+bz*dzvz
            Stretch= bx*Sx+by*Sy+bz*Sz
            bins_1 = np.geomspace( 1,1e9,19)
            bins_m1 = -bins_1[::-1]
            bins = nar(list(bins_m1)+list(bins_1))
            bincen = 0.5*(bins[1:]+bins[:-1])


            def splat(array, ax, title):
                smooth= ndimage.gaussian_filter1d(array, 2, 0)
                smooth = 0.5*(smooth[1:,:]+smooth[:-1,:])
                ds_x,ds_y,ds_h,ds_dv,ds_p=heat_map.heat_map( smooth.transpose(), tcenter, bins=bins, ax=ax)
                ax.set_yscale('symlog',linthresh=100)
                ax.set_title(title)
                return ds_x,ds_y,ds_h,ds_dv,ds_p

            splat(dxvx, ax9, 'dxvx')
            splat(dxvy, ax10, 'dxvy')
            splat(dxvz, ax11, 'dxvz')

            splat(dyvx, ax12,  'dyvx')
            splat(dyvy, ax13, 'dyvy')
            splat(dyvz, ax14, 'dyvz')

            splat(dzvx, ax15,  'dyvx')
            splat(dzvy, ax16, 'dyvy')
            splat(dzvz, ax17, 'dyvz')
            smooth_b=ndimage.gaussian_filter1d(B, 2, 0)
            db_dt = (smooth_b[1:,:]-smooth_b[:-1,:])/dt_square
            db_x,db_y,db_h,db_dv,db_p=heat_map.heat_map( db_dt.transpose(), tcenter, bins=bins, ax=ax0)
            ax0.set_yscale('symlog',linthresh=100)
            ax0.set_title('db/dt')

            Pdb1=splat( -B2*divv, ax2, '-B^2 divV')
            splat( Sx, ax1, 'Sx')
            splat( Sy, ax3, 'Sy')
            splat( Sz, ax5, 'Sz')
            Pdb2=splat( Stretch, ax4, 'Stretch')

            fig7,aaa=plt.subplots(2,2)
            aaax7=aaa[0][1]; aaa8=aaa[0][0]; aaa9=aaa[1][0];aaa10=aaa[1][1]
            aaax7.plot( bincen, Pdb2[2][0,:])
            aaax7.hist( Stretch[0,:], bins=bins)
            aaax7.set_xscale('symlog',linthresh=100)

            aaa8.plot( sorted( np.abs(Stretch[0,:])))
            aaa8.set(xlabel='nparticle', ylabel='Stretch')
            aaa8.set_yscale('log')

            aaa9.plot( sorted( np.abs(Stretch[0,:])), np.linspace(0,1,Stretch[0,:].size), marker='*')
            aaa9.set(ylabel='cuml',xlabel='||Stretch||', xscale='log')

            aaa10.hist( np.abs(Stretch[0,:]), bins=np.geomspace(100,1e9,19))
            print(np.abs(Stretch[0,:]).min())
            aaa10.set(xscale='log')
            #aaa10.set_xscale('symlog',linthresh=100)


            fig7.savefig('plots_to_sort/derp.png')

            def not_dumb_ticks(qaxis, labs, nticks, format="%0.1f"):
                ticks = qaxis.get_ticklocs()
                points = np.linspace( ticks.min(), ticks.max(), nticks)
                fake_points = np.linspace( ticks.min(), ticks.max(), labs.size)
                new_ticks = np.interp( points, fake_points, labs)
                qaxis.set_ticks(points)
                qaxis.set_ticklabels([format%num for num in new_ticks])


            fig.savefig('plots_to_sort/b_and_rho_3_c%04d.pdf'%(core_id))
            fig2,aaax=plt.subplots(1,3)
            moo = np.log10(db_h.transpose())
            dumb = np.zeros_like(moo)
            mork = np.log10(Pdb1[2].transpose())
            pork = np.log10(Pdb2[2].transpose())
            oot=np.stack([moo,mork,pork],axis=2)
            #oot=np.stack([moo,mork+pork,mork+pork],axis=2)
            aaax[0].imshow(oot, origin='lower', extent=[0,1,0,1])



                
            not_dumb_ticks(aaax[0].xaxis, times, 5, format="%0.3f")
            not_dumb_ticks(aaax[0].yaxis, nar([-1e9,-1e4,0,1e4,1e9]), 5, format="%0.1e")
            #aaax[0].xaxis.set_ticklabels(['arse','barse','carse'])
            aaax[0].set_title('dB, db1')

            oot=np.stack([moo,mork+pork, mork+pork],axis=2)
            #oot=np.stack([moo,mork+pork,mork+pork],axis=2)
            aaax[1].imshow(oot, origin='lower', extent=[0,1,0,1])

            venn1 = np.zeros([800,800])
            venn2 = np.zeros([800,800])
            venn3 = np.zeros([800,800])
            x,y = np.mgrid[0:1:1/800, 0:1:1/800]-0.5
            r=0.25
            x1,y1 = r*np.cos(2*np.pi/3), r*np.sin(2*np.pi/3)
            ok = (x-x1)**2 + (y-y1)**2 < 1.2*r**2
            
            venn1[ok] = 1
            x2,y2 = r*np.cos(4*np.pi/3), r*np.sin(4*np.pi/3)
            ok = (x-x2)**2 + (y-y2)**2 < 1.2*r**2
            venn2[ok] = 1
            x3,y3 = r*np.cos(6*np.pi/3), r*np.sin(6*np.pi/3)
            ok = (x-x3)**2 + (y-y3)**2 < 1.2*r**2
            venn3[ok] = 1
            venn = np.stack([venn1,venn2,venn3],axis=2)
            aaax[2].imshow(venn,origin='lower', extent=[-0.5,0.5,-0.5,0.5])
            aaax[2].text( y1,x1, 'dBdt', color='white')
            aaax[2].text( y2,x2, 'B2', color='white')
            aaax[2].text( y3,x3, 'Stretch', color='white')
            fig2.savefig('plots_to_sort/db_dt_c%04d.png'%core_id)
            return Pdb1

sim_list=['u902']
for sim in sim_list:
    ddd = dq_dt(TL.loops[sim])
    core_list=None
    #core_list=[7]
    core_list=[74]
    P=ddd.run(core_list=core_list)


