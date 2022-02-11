
from starter2 import *


sim = 'u16'; directory='/scratch1/dcollins/Paper19/u16b_sphere2_isolated'
#sim = 'u23'; directory='/scratch1/dcollins/Paper19/u23_point'
sim = 'u24'; directory='/scratch1/dcollins/Paper19/u24_point_periodic'
setname = "%s/DD%04d/data%04d"

frame=1
if 'ds' not in dir():
    ds = yt.load(setname%(directory,frame,frame))
if 1:
    yt.ProjectionPlot(ds,0,'density').save('plots_to_sort/%s_n%04d'%(sim,frame))
    yt.ProjectionPlot(ds,0,'GravPotential').save('plots_to_sort/%s_n%04d'%(sim,frame))

    ad = ds.all_data()
    print(ad['cell_mass'].sum())

    M = (ad['cell_mass'][ ad['density']>2]).sum()

foot = (0.5,0.5)
ray = ds.ortho_ray(0,foot)
plt.clf()
xc,yc,zc=[0.5]*3
r = np.sqrt((ray['x'].v-xc)**2 + (ray['y'].v-yc)**2 + (ray['z'].v-zc)**2)
G = ds['GravitationalConstant']/(np.pi*4)
rsphere=1./32
rsphere=0.25
#M = 4./3*np.pi*10*rsphere**3
if sim == 'u16':
    PhiA = -G*M*(3*rsphere**2-r**2)/(2*rsphere**3)
elif sim == 'u23':
    PhiA = np.zeros_like( r)
    ok1 = r>rsphere
    PhiA[ok1] = -G*M/r[ ok1]
    PhiA[~ok1] = -G*M*(3*rsphere**2-r[~ok1]**2)/(2*rsphere**3)

PhiA = np.zeros_like( r)
ok1 = r>rsphere
PhiA[ok1] = -G*M/r[ ok1]
PhiA[~ok1] = -G*M*(3*rsphere**2-r[~ok1]**2)/(2*rsphere**3)

#PhiA*=1.5
def dummy(x,m,b):
    if hasattr(PhiA,'v'):
        Phi = PhiA.v
    else:
        Phi = PhiA
    return Phi*m+b

from scipy.optimize import curve_fit



Phi1=ray['GravPotential']

popt,pcov = curve_fit( dummy, r, Phi1.v)
plt.plot( r, Phi1, label='data')
#rat=PhiA/Phi1
#plt.plot( r, rat)
rmin = np.argmin(r)
fac=1.5#Phi1[rmin]/PhiA[rmin]
plt.plot( r, PhiA, label='Analytic')
plt.plot( r, dummy( r, *popt), label='shift %0.1f %0.1f'%(popt[0],popt[1]))
plt.legend(loc=0)
plt.savefig('plots_to_sort/ray.png')
