from starter2 import *

from scipy import ndimage
import progressbar
import time
import movie_frames
reload(movie_frames)
import pcolormesh_helper as pch

class dq_dt2():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.collector=defaultdict(list)
    def run(self,core_list=None,frame_list=None, OOM=False):
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
            self.tcenter = 0.5*(times[:-1]+times[1:])
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
            self.collector['cores_used'].append(core_id)

            if True: #ms.nparticles < 1000:
                sl=slice(None)
                c=[0.5]*4
            else:
                sl = slice(None,None,10)
                #c=[0,0,0,0.1]
                c=[0.1]*4
            rho = ms.density[sl].transpose()
            dv = ms.cell_volume[sl].transpose()[mask,:]
            rho = rho[mask,:]
            Bmag=thtr.c([core_id],'magnetic_field_strength')[sl].transpose()[mask,:]#/colors.mean_field[this_looper.sim_name]
            BP = Bmag**2/2
            B2 = Bmag**2
            divv=thtr.c([core_id],'velocity_divergence')[sl].transpose()[mask,:]#/colors.mean_field[this_looper.sim_name]
            self.divv=divv


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

            if OOM:
                self.dxvx=dxvx
                self.dxvy=dxvy
                self.dxvz=dxvz
                self.dyvx=dyvx
                self.dyvy=dyvy
                self.dyvz=dyvz
                self.dzvx=dzvx
                self.dzvy=dzvy
                self.dzvz=dzvz
                self.dixj=[[self.dxvx,self.dxvy,self.dxvz],
                           [self.dyvx,self.dyvy,self.dyvz],
                           [self.dzvx,self.dzvy,self.dzvz]]
            Sx = bx*dxvx+by*dyvx+bz*dzvx
            Sy = bx*dxvy+by*dyvy+bz*dzvy
            Sz = bx*dxvz+by*dyvz+bz*dzvz
            Stretch= bx*Sx+by*Sy+bz*Sz

            R = Stretch/(B2*divv)
            self.R=R

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
            self.A0 = np.zeros_like(R)
            self.A1 = np.zeros_like(R)
            self.A2 = np.zeros_like(R)
            self.Br0 = np.zeros_like(R)
            self.Br1 = np.zeros_like(R)
            self.Br2 = np.zeros_like(R)
            self.Bp0 = np.zeros_like(R)
            self.Bp1 = np.zeros_like(R)
            self.Bp2 = np.zeros_like(R)

            self.theta=np.zeros_like(R)

            self.rho = rho
            self.dv = dv
            self.bmag = Bmag
            self.B2=B2

            self.B1dotV1 = np.zeros_like(R)
            self.B2dotV1 = np.zeros_like(R)
            self.dumb_test=np.zeros_like(R)
            bar = progressbar.ProgressBar(maxval=len(frames))
            bar.start()


            angle_dot=[]
            angle_guess=[]
            theta_1=[]
            theta_5=[]
            theta_4=[]
            junk=defaultdict(list)
            for nf,frame in enumerate(frames):
                #print('frame',frame)
                bar.update(nf)


                #nf = np.where(frame == self.this_looper.tr.frames)[0][0]

                for ip,ppp in enumerate(particles):
                    Dij = np.array([[dxvx[nf,ip], dxvy[nf,ip], dxvz[nf,ip]] ,
                                    [dyvx[nf,ip], dyvy[nf,ip], dyvz[nf,ip]] ,
                                    [dzvx[nf,ip], dzvy[nf,ip], dzvz[nf,ip]] ])/divv[nf,ip]
                    Bi =  np.array([bx[nf,ip],by[nf,ip],bz[nf,ip]])/Bmag[nf,ip]

                    U,H = scipy.linalg.polar(Dij)
                    A,E = scipy.linalg.eig(H)
                    #print(A.sum())
                    #print(Dij[0][0]+Dij[1][1]+Dij[2][2])
                    #print(np.trace(U))
                    #print(np.linalg.det(Dij))
                    #print(np.linalg.det(U))
                    if (A<=0).any() :
                        print("WTF negative eigen")
                    Lambda=np.diag(A)

                    b_rot = np.matmul(Bi,U)
                    b_rot_new = b_rot@E
                    b_new = E.transpose()@Bi
                    s_new = b_rot_new@Lambda
                    r_new_basis=np.dot(b_new,s_new)
                    #print(R[nf,ip]/r_new_basis)
                    Isrt = np.argsort(A)[::-1]
                    self.A0[nf,ip]=A[Isrt[0]].real
                    self.A1[nf,ip]=A[Isrt[1]].real
                    self.A2[nf,ip]=A[Isrt[2]].real
                    self.Br0[nf,ip]=b_new[Isrt[0]]
                    self.Br1[nf,ip]=b_new[Isrt[1]]
                    self.Br2[nf,ip]=b_new[Isrt[2]]

                    self.Bp0[nf,ip]=b_rot_new[Isrt[0]] 
                    self.Bp1[nf,ip]=b_rot_new[Isrt[1]]
                    self.Bp2[nf,ip]=b_rot_new[Isrt[2]]

                    self.theta[nf,ip]=(np.trace(U)-1)/2

                    if 1:
                        VVV = Bi
                        if 0:
                            #hijack U
                            from numpy import cross, eye, dot
                            from scipy.linalg import expm, norm
                            def M(axis, theta):
                                return expm(cross(eye(3), axis/norm(axis)*theta))

                            v, axis, theta = nar([1,0,0]), nar([0,1,0]), 0.8
                            #v=Bi
                            #v=v/(v*v).sum()**0.5
                            U = M(axis, theta)

                            #math I got to work the first time.
                            #axis=Arot
                            #ahat = axis/(axis*axis).sum()
                            #v1 = (ahat*v).sum()*ahat
                            #v2 = v-v1
                            #v2mag = np.sqrt((v2*v2).sum())
                            #Dprime = v2mag*np.sin( theta/2)
                            #vmag = np.sqrt((v*v).sum())
                            #theta2 = 2*np.arcsin( Dprime/vmag)
                            #print("is 0.8",theta2, 'isnt',Theta2, 'trace',Theta1)


                        au,eu=np.linalg.eig(U)
                        primary_axis = np.argmax(np.abs(au.real)) #the one thats 1.
                        Arot = eu[:,primary_axis]+0
                        Arot /= (Arot*Arot).sum()**0.5
                        if (np.sum(Arot.imag)) /(np.sum(Arot.real)) > 1e-7:
                            print("imaginary eigen vector")
                            pdb.set_trace()
                        Arot = Arot.real
                        #print(eu)
                        #print(au)
                        #print(primary_axis)
                        #print(np.abs(au.real))


                        V = VVV #Bi+0
                        Vm = (V*V).sum()**0.5
                        V1 = (Arot*V).sum()*Arot
                        V2 = V-V1
                        V2m = (V2*V2).sum()**0.5
                        V1m = (V1*V1).sum()**0.5
                        if 0:
                            print('kludge')
                            V = V2
                            Vm = (V*V).sum()**0.5

                        #the other one, not necessary
                        #next_axis = (primary_axis+2)%3
                        #next_eig = au[next_axis]
                        #theta_5.append( np.angle(next_eig))

                        next_axis = (primary_axis+1)%3
                        next_eig = au[next_axis]
                        Theta4 = np.angle(next_eig)
                        Theta1 = np.arccos((np.trace(U)-1)/2)
                        #Theta2 = 2*np.arcsin( V2m/Vm*np.sin(Theta4/2))
                        if au[primary_axis]>0:
                            Theta2 = 2*np.arcsin( V2m/Vm*np.sin(Theta4/2))
                        else:
                            thisthing= V2m*V2m*np.cos(Theta4)-V1m**2
                            Theta2 = np.arccos(thisthing)

                        theta_1.append(Theta1)
                        theta_4.append(Theta4)

                        #V=V2
                        VR = U@V
                        VRm = (VR*VR).sum()**0.5
                        Theta3 = np.arccos((V*VR).sum()/(VRm*Vm))

                        #print("A",Arot)
                        #print('d',(V2*Arot).sum())
                        #print("Trace %0.4f Guess %0.4f Dot %0.4f rat %0.4f"%(Theta1, Theta2, Theta3,Theta2/Theta3))


                        #if au[primary_axis]<0:
                        angle_guess.append(Theta2)
                        angle_dot.append(Theta3)
                        #    junk['t1'].append(Theta1)
                        #    junk['t4'].append(Theta4)
                        #    angle_guess.append(Theta2)
                        #    angle_dot.append(Theta3)
                        #    junk['neg'].append(1)
                        #else:
                        #    junk['pos'].append(1)


                        #Bi2=np.array([1,0,0])
                        #b2 = np.matmul(U,Bi2)
                        #Theta1 = (np.trace(U)-1)/2
                        ##Theta2 = b2[0]*Bi2[0]+b2[1]*Bi2[1]+b2[2]*Bi2[2]
                        #Theta2 = np.dot(b2,Bi2)/np.sqrt(np.dot(b2,b2)*np.dot(Bi2,Bi2))
                        #Theta3=np.cos(np.angle(au[1]))
                        ##Theta3=0.5*(au[1]+au[2])
                        ##Theta1 = Theta2+Theta3
                        #print("dU %0.1f T 1 %0.2f T2 %0.2f T3 %0.2f"%(np.linalg.det(U),Theta1,Theta2,Theta3))
                        if 'num' not in dir():
                            num=0
                        #if num > 10:
                        #    break
                        num+=1  

                    #iii0,iii1,iii2=0,1,2
                    iii0,iii1,iii2=Isrt
                    #print(Isrt)
                    r_eig=b_new[iii0]*b_rot_new[iii0]*A[iii0]+\
                          b_new[iii1]*b_rot_new[iii1]*A[iii1]+\
                          b_new[iii2]*b_rot_new[iii2]*A[iii2]
                    #print(R[nf,ip]/r_eig)
                    #self.arf.append(R[nf,ip]/r_eig)
                    #self.arf.append(r_eig)
                    self.dumb_test[nf,ip]=r_eig.real



                    E1 = E[:,0]
                    self.B1dotV1[nf,ip] = (E1*b_rot).sum()
                    self.B2dotV1[nf,ip] = (E1*b_rot_new).sum()


            #rho_0 = (self.dv*self.rho)[0,:].sum()/self.dv[0,:].sum()
            #B_0 = (self.dv*self.bmag)[0,:].sum()/self.dv[0,:].sum()
            rho_0=rho[0,:].mean()
            B_0=B2[0,:].mean()
            pfit = np.polyfit( np.log10(self.B2/B_0).flatten(), np.log10(self.rho/rho_0).flatten(),1)
            self.collector['pfit0'].append(pfit[0])
            self.collector['pfit1'].append(pfit[1])
            self.collector['rho_0'].append(rho_0)
            self.collector['B_0'].append(B_0)
            RB_Mean = (rho*R*dv).sum()/(rho*dv).sum()
            self.collector['RB_mean'].append(RB_Mean)

            self.core_id=core_id

            if 0:
                fig3,ax3=plt.subplots(2,2)
                ax3[0][0].scatter(angle_dot,angle_guess)
                ax3[0][1].scatter(theta_1,theta_4)
                #ax3[1][0].scatter(nar(junk['t1']),nar(junk['t4']))
                #pdb.set_trace()
                #ax3[1][0].set(yscale='log')
                x=nar(angle_dot)
                y=nar(angle_guess)
                pch.simple_phase(x,y,ax=ax3[1][1])
                fig3.savefig('plots_to_sort/angle_%s_c%04d'%(this_looper.sim_name,core_id))
