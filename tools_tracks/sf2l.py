
from starter2 import *
import scipy.signal
import matplotlib.patches as patches
plt.close('all')
figsize = None #(12,12)
import tools.radial_binner as rb
reload(rb)
from scipy.optimize import curve_fit
import tools_spectra.spectra_tools as st
reload(st)

#a = np.arange(27)
#a.shape=(3,3,3)
#print(a)
#print("=====")
#print(np.roll(a,(1,0,0),axis=(0,1,2)))
#print("=====")
#print(np.roll(a,(0,1,0),axis=(0,1,2)))
#print("=====")
#print(np.roll(a,(0,0,1),axis=(0,1,2)))
def get_sf2l(ad,subset=1):
    vx=ad['velocity_x']
    vy=ad['velocity_y']
    vz=ad['velocity_z']
    SF2 = np.zeros_like(vx.v)
    Rmag = np.zeros_like(SF2)
    Nx = SF2.shape[0]
    if subset == 0:
        Half = int(Nx)
        hunt_range=range(-Half,Half,4)
    elif subset == 1:
        Half = int(Nx//4)
        hunt_range=range(-Half,Half,4)
    elif subset == 2:
        Half = int(Nx//8)
        hunt_range=range(-Half,Half,4)
    else:
        print("Bad coding")
        raise

    for i in hunt_range:
        for j in hunt_range:
            for k in hunt_range:
                R = np.array([i,j,k])
                ThisRmag = np.sqrt((R*R).sum())
                if ThisRmag == 0:
                    continue
                Rhat = R/ThisRmag
                print(R, Rhat)
                dvx= vx - np.roll(vx,R,axis=(0,1,2))
                dvy= vy - np.roll(vy,R,axis=(0,1,2))
                dvz= vz - np.roll(vz,R,axis=(0,1,2))
                SF2[i,j,k]=( (dvx*Rhat[0]+dvy*Rhat[1]+dvz*Rhat[2])**2).sum()
                Rmag[i,j,k] = ThisRmag
    return Rmag,SF2


def make_rmag(ad):
    vx=ad['velocity_x']
    vy=ad['velocity_y']
    vz=ad['velocity_z']
    Rvec = np.zeros_like(vx)
    hunt_range=range(0,128,4)
    for i in hunt_range:
        for j in hunt_range:
            for k in hunt_range:
                R = np.array([i,j,k])
                Rmag = np.sqrt((R*R).sum())
                if Rmag == 0:
                    continue
                print(i,j,k)
                Rvec[i,j,k] = Rmag
    return Rvec

