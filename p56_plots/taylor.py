
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

if 'ds' not in dir():
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    ds = yt.load("%s/DD%04d/data%04d"%(directory,5,5))
    left=[0.0]*3
    resolution = ds['TopGridDimensions'] 
    cg=ds.covering_grid(0,left,resolution,num_ghost_zones=2)
if 1:
    #rho1 = cg['density'].v #-1 #[:40,:40,:40]
    #rho2=rho1
    vx = cg['velocity_x'].v
    vy = cg['velocity_y'].v
    vz = cg['velocity_z'].v
    omega2 = cg['vorticity_magnitude']**2
    total_volume = cg['cell_volume'].sum()
#   KE = (cg['kinetic_energy']*cg['cell_volume']).sum()
    KE = (cg['velocity_magnitude']**2*cg['cell_volume']).sum()
    Omega = 0.5*(omega2*cg['cell_volume']).sum()
    Taylor = np.sqrt(5*KE/Omega)
#   print(Taylor)

#   sigma_vx = np.std(vx)
#   sigma_vy = np.std(vy)
#   sigma_vz = np.std(vz)

