
from starter2 import *
import scipy.signal
import matplotlib.patches as patches
plt.close('all')
figsize = None #(12,12)
import tools.radial_binner as rb
reload(rb)
from scipy.optimize import curve_fit
def powerlaw(r,rho0, r0, alpha):
    return alpha*np.log10(r/r0) + np.log10(rho0)
def powerlaw2(r, r0,alpha):
    rhosquared=1
    return alpha*np.log10(r/r0) + np.log10(rhosquared)
axis=0
twopi = np.pi*2

import three_loopers_1tff as tl

Taylor_v = {} #velocity only
Taylor_k = {} #kinetic energy
Taylor_d = {} #with divergence.
for loop in [tl.looper1, tl.looper2, tl.looper3]:
    #directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    #ds = yt.load("%s/DD%04d/data%04d"%(directory,5,5))
    name = loop.out_prefix
    ds = loop.load(5)
    left=[0.0]*3
    resolution = ds['TopGridDimensions'] 
    cg=ds.covering_grid(0,left,resolution,num_ghost_zones=2)
    #rho1 = cg['density'].v #-1 #[:40,:40,:40]
    #rho2=rho1
    vx = cg['velocity_x'].v
    vy = cg['velocity_y'].v
    vz = cg['velocity_z'].v
    omega2 = cg['vorticity_magnitude']**2
    divv2 = cg['velocity_divergence']**2
    total_volume = cg['cell_volume'].sum()

    KE = (cg['velocity_magnitude']**2*cg['cell_volume']).sum()
    Omega = 0.5*(omega2*cg['cell_volume']).sum()
    Taylor_v[name] = np.sqrt(5*KE/Omega)
    print("Taylor v",Taylor_v[name]*128)

    KE2 = (cg['kinetic_energy']*cg['cell_volume']).sum()
    O2  = (cg['density']*omega2*cg['cell_volume']).sum()
    Taylor_k[name] = np.sqrt(5*KE2/O2)
    print("Traylor K",Taylor_k[name]*128)


    KE2 = (cg['kinetic_energy']*cg['cell_volume']).sum()
    O2  = (cg['density']*omega2*cg['cell_volume']).sum()
    D2  = (cg['density']*divv2*cg['cell_volume']).sum()
    Taylor_d[name] = np.sqrt(5*KE2/(O2+4./3*D2))
    print("Traylor D",Taylor_d[name]*128)

    
#   KE = (cg['kinetic_energy']*cg['cell_volume']).sum()

#   sigma_vx = np.std(vx)
#   sigma_vy = np.std(vy)
#   sigma_vz = np.std(vz)

