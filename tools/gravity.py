from starter2 import *
import xtra_energy

class gravity_complex():
    def __init__(self,density,G):
        self.density=density
        self.G=G

    def solve(self):
        rho = self.density.v
        rhohat_me = np.fft.fftn(rho)
        nxy = rhohat_me.shape
        #khat = np.mgrid[-nxy[0]/2:nxy[0]/2:1, -nxy[1]/2:nxy[1]/2:1, 0:nxy[2]:1]*np.pi*2
        #k2=khat[0,...]**2+khat[1,...]**2+khat[2,...]**2
        ka = np.fft.fftfreq(nxy[0])*nxy[0]*np.pi*2
        kb = np.fft.fftfreq(nxy[0])*nxy[0]*np.pi*2
        kc = np.fft.fftfreq(nxy[0])*nxy[0]*np.pi*2
        khat = np.meshgrid(ka,kb,kc,indexing='ij')
        k2=khat[0]**2+khat[1]**2+khat[2]**2
        #Real multi-d 'real' fft is complicated.
        #Last axis is a real transform, and only half as long.
        #The other two need to be rotated.

        #khat0 = np.mgrid[-nxy[0]/2:nxy[0]/2:1, -nxy[1]/2:nxy[1]/2:1, -nxy[2]/2:nxy[2]/2:1]*np.pi*2
        #k20=khat0[0,...]**2+khat0[1,...]**2+khat0[2,...]**2
        #k20 = np.fft.fftshift(k20,axes=0)
        #k20 = np.fft.fftshift(k20,axes=1)
        #k20 = np.fft.fftshift(k20,axes=2)
        k2inv = np.zeros_like(k2)
        k2inv[ k2>0] = 1/k2[k2>0]
        k2inv *= -self.G
        self.green = k2inv
        self.phihat = (rhohat_me*k2inv)
        self.phi = np.fft.ifftn(self.phihat)

        #store
        self.khat = khat
        self.k2inv=k2inv
        self.k2=k2
        self.rhohat=rhohat_me
        
        also_rhohat = self.phihat/k2inv
        also_rhohat[0,0,0]=self.density.sum().v
        self.also_rho = np.fft.ifftn(also_rhohat)

    def get_g(self):
        #complex
        khatx=self.khat[0]
        khaty=self.khat[1]
        khatz=self.khat[2]#real FFT, z direction is the only real one.
        self.khatx=khatx#kludge
        self.gx_hat = khatx*self.phihat
        self.gy_hat = khaty*self.phihat
        self.gz_hat = khatz*self.phihat
        self.dxgx_hat = khatx**2*self.phihat#*-1/self.G
        self.dygy_hat = khaty**2*self.phihat#*-1/self.G
        self.dzgz_hat = khatz**2*self.phihat#*-1/self.G
        total = self.density.sum()
        self.gx_hat[0,0,0]=total/3*-self.G
        self.gy_hat[0,0,0]=total/3*-self.G
        self.gz_hat[0,0,0]=total/3*-self.G
        self.dxgx_hat[0,0,0]=total/3*-self.G
        self.dygy_hat[0,0,0]=total/3*-self.G
        self.dzgz_hat[0,0,0]=total/3*-self.G
        print("x c",khatx[:2,:2,:2])
        print("y c",khaty[:2,:2,:2])
        self.gx = np.fft.ifftn(self.gx_hat)
        self.gy = np.fft.ifftn(self.gy_hat)
        self.gz = np.fft.ifftn(self.gz_hat)
        self.dxgx = np.fft.ifftn(self.dxgx_hat)
        self.dygy = np.fft.ifftn(self.dygy_hat)
        self.dzgz = np.fft.ifftn(self.dzgz_hat)

        #Ix = (np.abs(self.gx.imag)/np.abs(self.gx.real)).sum()
        #print("IMAGINARY",Ix)

        also_rhohat2 = (self.dxgx_hat+self.dygy_hat+self.dzgz_hat)
        #also_rhohat2 = (self.dxgx_hat+self.dygy_hat+self.dzgz_hat)
        #also_rhohat2[0,0,0]=self.density.sum().v
        #print("dood",self.also_also_rho.shape)
        self.also_also_rho=np.fft.ifftn(also_rhohat2)
        self.also_also_also_rho = (self.dxgx+self.dygy+self.dzgz)/-self.G
if 0:
        also_rhohat = self.phihat/self.k2inv
        also_rhohat[0,0,0]=also_rhohat.size
        #print("RHO 0", rhohat_me[0,0,0])
        print(np.abs(also_rhohat2-also_rhohat).sum())
        self.derp=np.fft.irfftn(also_rhohat)
        #print("rho error", np.abs(rhohat_me-also_rhohat).sum()/also_rhohat.size)
        self.also_rho = np.fft.irfftn(also_rhohat)
        #self.also_rho = np.fft.irfftn(rhohat_me)
        #self.also_rho = rho
        #print("k rho error" , np.abs(rho - self.also_rho).sum())


class gravity_real():
    def __init__(self,density,G):
        self.density=density
        self.G=G

    def solve(self):
        rho = self.density.v
        rhohat_me = np.fft.rfftn(rho)
        self.rhohat=rhohat_me
        nxy = rhohat_me.shape
        khat = np.mgrid[-nxy[0]/2:nxy[0]/2:1, -nxy[1]/2:nxy[1]/2:1, 0:nxy[2]:1]*np.pi*2
        k2=khat[0,...]**2+khat[1,...]**2+khat[2,...]**2
        #Real multi-d 'real' fft is complicated.
        #Last axis is a real transform, and only half as long.
        #The other two need to be rotated.
        k2 = np.fft.fftshift(k2,axes=0)
        k2 = np.fft.fftshift(k2,axes=1)
        k2inv = np.zeros_like(k2)
        k2inv[ k2>0] = 1/k2[k2>0]
        k2inv *= -self.G
        self.green = k2inv
        self.phihat = (rhohat_me*k2inv)
        self.phi = np.fft.irfftn(self.phihat)

        #store
        self.khat = khat
        self.k2inv=k2inv
        self.k2=k2
        
        also_rhohat = self.phihat/k2inv
        also_rhohat[0,0,0]=self.density.sum().v
        self.also_rho = np.fft.irfftn(also_rhohat)

    def get_g(self):
        #real
        khatx=self.khat[0]
        khatx=np.fft.fftshift(khatx,axes=0)
        khatx=np.fft.fftshift(khatx,axes=1)
        khaty=self.khat[1]
        khaty=np.fft.fftshift(khaty,axes=0)
        khaty=np.fft.fftshift(khaty,axes=1)
        khatz=self.khat[2]#real FFT, z direction is the only real one.
        self.khatx=khatx#kludge
        #khatz=np.fft.fftshift(khatz,axes=0)
        #khatz=np.fft.fftshift(khatz,axes=1)
        self.gx_hat = khatx*self.phihat
        self.gy_hat = khaty*self.phihat
        self.gz_hat = khatz*self.phihat
        self.dxgx_hat = khatx**2*self.phihat#*-1/self.G
        self.dygy_hat = khaty**2*self.phihat#*-1/self.G
        self.dzgz_hat = khatz**2*self.phihat#*-1/self.G
        print("x r",khatx[:2,:2,:2])
        print("y r",khaty[:2,:2,:2])
        total = self.density.sum()
        self.gx_hat[0,0,0]=total/3*-self.G
        self.gy_hat[0,0,0]=total/3*-self.G
        self.gz_hat[0,0,0]=total/3*-self.G
        self.dxgx_hat[0,0,0]=total/3*-self.G
        self.dygy_hat[0,0,0]=total/3*-self.G
        self.dzgz_hat[0,0,0]=total/3*-self.G
        self.gx = np.fft.irfftn(self.gx_hat)
        self.gy = np.fft.irfftn(self.gy_hat)
        self.gz = np.fft.irfftn(self.gz_hat)
        self.dxgx = np.fft.irfftn(self.dxgx_hat)
        self.dygy = np.fft.irfftn(self.dygy_hat)
        self.dzgz = np.fft.irfftn(self.dzgz_hat)

        also_rhohat2 = (self.dxgx_hat+self.dygy_hat+self.dzgz_hat)
        #also_rhohat2 = (self.dxgx_hat+self.dygy_hat+self.dzgz_hat)
        #also_rhohat2[0,0,0]=self.density.sum().v
        #print("dood",self.also_also_rho.shape)
        self.also_also_rho=np.fft.irfftn(also_rhohat2)
        self.also_also_also_rho = (self.dxgx+self.dygy+self.dzgz)/-self.G
if 0:
        also_rhohat = self.phihat/self.k2inv
        also_rhohat[0,0,0]=also_rhohat.size
        #print("RHO 0", rhohat_me[0,0,0])
        print(np.abs(also_rhohat2-also_rhohat).sum())
        self.derp=np.fft.irfftn(also_rhohat)
        #print("rho error", np.abs(rhohat_me-also_rhohat).sum()/also_rhohat.size)
        self.also_rho = np.fft.irfftn(also_rhohat)
        #self.also_rho = np.fft.irfftn(rhohat_me)
        #self.also_rho = rho
        #print("k rho error" , np.abs(rho - self.also_rho).sum())
