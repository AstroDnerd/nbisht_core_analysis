from starter2 import *

import three_loopers_u900 as TL

from scipy import ndimage

import movie_frames
reload(movie_frames)
import heat_map
reload(heat_map)
plt.close('all')
import pcolormesh_helper as pch
reload(pch)

def splat(array, tcenter, ax, title,bins):
    smooth= ndimage.gaussian_filter1d(array, 2, 0)
    smooth = 0.5*(smooth[1:,:]+smooth[:-1,:])
    ds_x,ds_y,ds_h,ds_dv,ds_p=heat_map.heat_map( smooth.transpose(), tcenter, bins=bins, ax=ax)
    ax.set_yscale('symlog',linthresh=100)
    ax.set_title(title)
    return ds_x,ds_y,ds_h,ds_dv,ds_p
class dq_dt2():
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
            dt2 = times[2:]-times[:-2]
            dt2.shape = dt2.size,1
            tcen2 = 0.5*(times[2:]+times[:-2])

            dt_square = dt+0
            dt_square.shape = dt_square.size,1


        if core_list is None:
            core_list = sorted(np.unique( thtr.core_ids))


        for core_id in core_list:
            print('go',core_id)
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            ms.particle_pos(core_id)

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
            BP = Bmag**2/2
            B2 = Bmag**2
            divv=thtr.c([core_id],'velocity_divergence')[sl].transpose()[mask,:]#/colors.mean_field[this_looper.sim_name]

            fig, ax=plt.subplots(3,4,figsize=(20,12))
            ax0=ax[0][0]
            ax1=ax[0][1]
            ax2=ax[0][2]
            ax3=ax[1][0]
            ax4=ax[1][1]
            ax5=ax[1][2]
            ax6=ax[2][0]
            ax7=ax[2][1]
            ax8=ax[2][2]
            if 0:
                ax0=ax[0][0];
                ax1=ax[0][1] 
                ax2=ax[0][2];
                #ax3=ax[0][3]
                #ax4=ax[0][4];
                #ax5=ax[0][5]

                ax6=ax[1][0];
                ax7=ax[1][1]
                ax8=ax[1][2]
                ax9=ax[1][3]
                ax10=ax[1][4]
                ax11=ax[1][5]

                ax12=ax[2][0]
                ax13=ax[2][1]
                ax14=ax[2][2]
                ax15=ax[2][3]
                ax16=ax[2][4]
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

            R = Stretch/(B2*divv)

            SxN = Sx/np.sqrt(B2)/divv
            SyN = Sy/np.sqrt(B2)/divv
            SzN = Sz/np.sqrt(B2)/divv

            bins_1 = np.geomspace( 1,1e9,19)
            bins_m1 = -bins_1[::-1]
            bins = nar(list(bins_m1)+list(bins_1))
            bincen = 0.5*(bins[1:]+bins[:-1])

            frames=np.arange( R.shape[0])
            particles=np.arange( R.shape[1])
            #particles=np.arange(R.shape[1]/4) #kludge
            A0 = np.zeros_like(R)
            A1 = np.zeros_like(R)
            A2 = np.zeros_like(R)
            Br0 = np.zeros_like(R)
            Br1 = np.zeros_like(R)
            Br2 = np.zeros_like(R)
            Bp0 = np.zeros_like(R)
            Bp1 = np.zeros_like(R)
            Bp2 = np.zeros_like(R)

            B1dotV1 = np.zeros_like(R)
            B2dotV1 = np.zeros_like(R)

            collector=[]
            for nf,frame in enumerate(frames):
                print('frame',frame)

                #nf = np.where(frame == self.this_looper.tr.frames)[0][0]

                for ip,ppp in enumerate(particles):
                    Dij = np.array([[dxvx[nf,ip], dxvy[nf,ip], dxvz[nf,ip]] ,
                                    [dyvx[nf,ip], dyvy[nf,ip], dyvz[nf,ip]] ,
                                    [dzvx[nf,ip], dzvy[nf,ip], dzvz[nf,ip]] ])/divv[nf,ip]
                    Bi =  np.array([bx[nf,ip],by[nf,ip],bz[nf,ip]])/Bmag[nf,ip]

                    U,H = scipy.linalg.polar(Dij)
                    A,E = scipy.linalg.eig(H)
                    if (A<=0).any() :
                        print("WTF negative eigen")
                    Lambda=np.diag(A)

                    b_rot = np.matmul(Bi,U)
                    b_rot_new = b_rot@E
                    b_new = E.transpose()@Bi
                    s_new = b_rot_new@Lambda
                    r_new_basis=np.dot(b_new,s_new)

                    Isrt = np.argsort(A)[::-1]
                    A0[nf,ip]=A[Isrt[0]]
                    A1[nf,ip]=A[Isrt[1]]/A[Isrt[0]]
                    A2[nf,ip]=A[Isrt[2]]/A[Isrt[0]]
                    Br0[nf,ip]=b_rot[Isrt[0]]
                    Br1[nf,ip]=b_rot[Isrt[1]]
                    Br2[nf,ip]=b_rot[Isrt[2]]

                    Bp0[nf,ip]=b_rot_new[Isrt[0]] #kludge
                    Bp1[nf,ip]=b_rot_new[Isrt[1]]
                    Bp2[nf,ip]=b_rot_new[Isrt[2]]


                    E1 = E[:,0]
                    B1dotV1[nf,ip] = (E1*b_rot).sum()
                    B2dotV1[nf,ip] = (E1*b_rot_new).sum()
                    #print('----')
                    #kprint( np.abs((H@E1/A[0]).real-E1).sum())




                    #print( "Brot", (b_rot**2).sum())
                    #print( "BrotN", (b_rot_new**2).sum())
                    #collector.append( A.sum())
                    #collector.append(R[nf,ip]/r_new_basis)
            #collector=nar(collector)
            #return collector
            #print(collector.imag)
            #print(collector.real.mean(), (collector.real-collector.real.mean()).sum())
            #pdb.set_trace()
            ext=extents()
            #ext(A0);ext(A1);ext(A2)
            ext(A0[A0>0]);ext(A1[A1>0]);ext(A2[A2>0])
            #print(A0.min())
            bins2=np.geomspace(ext.minmax[0],ext.minmax[1],64)

            splat(A0, tcenter, ax0, 'A0',bins=bins2)
            splat(A1, tcenter, ax1, 'A1/A0',bins=bins2)
            splat(A2, tcenter, ax2, 'A2/A0',bins=bins2)
            ax0.set(yscale='log',ylim=ext.minmax)
            ax1.set(yscale='log',ylim=ext.minmax)
            ax2.set(yscale='log',ylim=ext.minmax)

            def sp2(arr,t,ax,lab):
                b2 = np.abs(arr)
                bins=np.geomspace(b2[b2>0].min(),b2.max(),64)
                smooth=arr+0
                smooth = 0.5*(smooth[1:,:]+smooth[:-1,:])
                ds_x,ds_y,ds_h,ds_dv,ds_p=heat_map.heat_map( smooth.transpose(), tcenter, bins=bins, ax=ax)
                ax.set_yscale('log')
            def sp3(arr,t,ax,lab):
                for n in np.arange( arr.shape[0]):
                    b=arr[n,:]
                    b.sort()
                    print(np.abs(b).min())
                    #b=np.log10(np.abs(b))

                    y = np.arange(b.size)/b.size
                    #ax.plot(b,y)
                    ax.plot(b)

            #sp3(Br0,tcenter,ax3,'Bp0')
            #sp2(Bp0,tcenter,ax3,'Bp0')
            if 0:
                #B vs B phases.  Makes circles.
                #pch.simple_phase( np.log10((Bp0**2).flatten()), np.log10((Bp1**2).flatten()), ax=ax3)
                pch.simple_phase( ((Bp0).flatten()), ((Bp1).flatten()), ax=ax3)
                ax3.set_aspect('equal')
                pch.simple_phase( Bp0.flatten(), Bp2.flatten(), ax=ax4)
                ax4.set_aspect('equal')
                pch.simple_phase( ((Br0).flatten()), ((Br1).flatten()), ax=ax6)
                ax6.set_aspect('equal')
                pch.simple_phase( Br0.flatten(), Br2.flatten(), ax=ax7)
                ax7.set_aspect('equal')
            if 1:

                ext=extents()
                ext(Br0);ext(Br1);ext(Br2)
                bins=np.linspace(ext.minmax[0],ext.minmax[1],64)

                splat(Bp0, tcenter, ax3, 'B_r_n 0',bins=bins)
                splat(Bp1, tcenter, ax4, 'B_r_n 1',bins=bins)
                splat(Bp2, tcenter, ax5, 'B_r_n 2',bins=bins)
                #ax3.set(yscale='linear',ylim=ext.minmax)
                #ax3.set_yscale('symlog',linthresh=1e-2)
                ax3.set(yscale='linear',ylim=ext.minmax)
                ax4.set(yscale='linear',ylim=ext.minmax)
                ax5.set(yscale='linear',ylim=ext.minmax)

                splat(Br0, tcenter, ax6, 'B_n 0',bins=bins)
                splat(Br1, tcenter, ax7, 'B_n 1',bins=bins)
                splat(Br2, tcenter, ax8, 'B_n 2',bins=bins)
                ax6.set(yscale='linear',ylim=ext.minmax)
                ax7.set(yscale='linear',ylim=ext.minmax)
                ax8.set(yscale='linear',ylim=ext.minmax)
            fig.savefig('plots_to_sort/the_stuff_%s_c%04d'%(self.this_looper.sim_name,core_id))
            if 0:

                    #print(R[nf,ip]-r_new_basis)

                    if 0:
                        if 0:
                            Toy = np.array([[0,1,0],
                                            [0,0,0],
                                            [0,0,0]])
                            print(np.matmul(Toy,Bi))
                            print(Bi)
                            print('---')

                        s_hat = np.matmul(Bi,Dij)
                        #print('Shat',s_hat)
                        DivV=divv[nf,ip]
                        print('SSSS',SxN[nf,ip],SyN[nf,ip],SzN[nf,ip])
                        print('RRRRRR',R[nf,ip])
                        #print('Bi',Bi)
                        #print('---')
                        if 0:
                            also_r = np.dot(Bi,s_hat)
                            print('Also R',also_r)
                            print('RRRRRR',R[nf,ip])
                            print('---')
                        U,H = scipy.linalg.polar(Dij)
                        A,E = scipy.linalg.eig(H)
                        Lambda=np.diag(A)
                        #print(A)

                        #works
                        b_rot = np.matmul(Bi,U)
                        S2 = np.matmul(b_rot,H)

                        print('---')
                        #Lambda=E.transpose()@H@E
                        #AlsoH =  E@Lambda@E.transpose()
                        #S_3= b_rot@E@Lambda@E.transpose()
                    




                    #print(S2)
                    #print(s_hat)



                    if 0:
                        #
                        # Just getting the eigen vectors doesn't work
                        #
                        print("--- Dij ---")
                        print(Dij)
                        print("Det Dij",scipy.linalg.det(Dij))
                        val,right=scipy.linalg.eig(Dij)
                        print("val",val, "sum",val.sum())
                        print("right",right)

                    if 0:
                        U,H = scipy.linalg.polar(Dij)
                        print("Det U",scipy.linalg.det(U))
                        print("Det H",scipy.linalg.det(H))
                        print("Val U", scipy.linalg.eig(U)[0])
                        print("Val H", scipy.linalg.eig(H)[0])

                    if 0:
                        #
                        # a aT is a projection operator.
                        # You can't rotate a projection operator, it only rotates the first vector
                        Bp = np.array([[ Bi[0]*Bi[0], Bi[0]*Bi[1],Bi[0]*Bi[2]],
                                       [ Bi[1]*Bi[0], Bi[1]*Bi[1],Bi[1]*Bi[2]],
                                       [ Bi[2]*Bi[0], Bi[2]*Bi[1],Bi[2]*Bi[2]]])

                        print("Bi",Bi)
                        Bj = np.matmul(U,Bi)


                        B2 = np.array([[ Bj[0]*Bi[0], Bj[0]*Bi[1],Bj[0]*Bi[2]],
                                       [ Bj[1]*Bi[0], Bj[1]*Bi[1],Bj[1]*Bi[2]],
                                       [ Bj[2]*Bi[0], Bj[2]*Bi[1],Bj[2]*Bi[2]]])

                        B3 = np.matmul(U,Bp)
                        print("B2",B2)
                        print("B3",B3)
                        print("B2-B3",B2-B3)







sim_list=['u902']
for sim in sim_list:
    ddd = dq_dt2(TL.loops[sim])
    core_list=None
    #core_list=[7]
    #core_list=[74]#, 112]
    core_list=[112]
    P=ddd.run(core_list=core_list)


