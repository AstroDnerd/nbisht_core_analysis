from starter2 import *
import xtra_energy


class gravity():
    def __init__(self,density,G):
        self.density=density
        self.G=G

    def solve(self):
        rho = self.density
        rhohat_me = np.fft.rfftn(rho)
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
        self.phihat = rhohat_me*k2inv
        self.phi = np.fft.irfftn(self.phihat)

