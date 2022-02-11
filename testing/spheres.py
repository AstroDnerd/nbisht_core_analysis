
from starter2 import *


if 0:
    sim = 'u16'; 
    directory='/scratch1/dcollins/Paper19/u16b_sphere2_isolated'
    rsphere=0.25
#sim = 'u23'; directory='/scratch1/dcollins/Paper19/u23_point'
if 0:
    sim = 'u24'; 
    directory='/scratch1/dcollins/Paper19/u24_point_periodic'
    rsphere=1./32
if 1:
    sim = 'u25'; 
    directory='/scratch1/dcollins/Paper19/u25_little_big'
    rsphere=0.25
setname = "%s/DD%04d/data%04d"

frame=1
if 'ds' not in dir():
    ds = yt.load(setname%(directory,frame,frame))

if 1:
    #make some projections
    yt.ProjectionPlot(ds,0,YT_density).save('plots_to_sort/%s_n%04d'%(sim,frame))
    yt.ProjectionPlot(ds,0,YT_potential).save('plots_to_sort/%s_n%04d'%(sim,frame))

    ad = ds.all_data()
    print(ad['cell_mass'].sum())

    M = (ad['cell_mass'][ ad['density']>2]).sum()

if 0:
    #a ray from the origin
    #Use either this or the next one
    foot = (0.5,0.5)
    ray = ds.ortho_ray(0,foot)
    plt.clf()
    xc,yc,zc=[0.5]*3
    r = np.sqrt((ray[YT_x].v-xc)**2 + (ray[YT_y].v-yc)**2 + (ray[YT_z].v-zc)**2)
    G = ds['GravitationalConstant']/(np.pi*4)
    Phi_data=ray['GravPotential']

if 1:
    #all the points in the sphere
    ad = ds.all_data()
    xc,yc,zc=[0.5]*3
    r = np.sqrt( (ad[YT_x].v-xc)**2+(ad[YT_y].v-yc)**2+(ad[YT_z].v-zc)**2)
    Phi_data = ad[YT_potential]


#analytic result
#Dan, make sure you understand what I'm doing here.
Phi_analytic = np.zeros_like( r)
ok1 = r>rsphere
Phi_analytic[ok1] = -G*M/r[ ok1]
Phi_analytic[~ok1] = -G*M*(3*rsphere**2-r[~ok1]**2)/(2*rsphere**3)


def shifter(x,m,b):
    #Shifter slides the analytic solution Phi_analytic around to match the data
    if hasattr(Phi_analytic,'v'):
        Phi = Phi_analytic.v
    else:
        Phi = Phi_analytic
    return Phi*m+b

from scipy.optimize import curve_fit

#fit the analytic, which is stored in the function shifter
popt,pcov = curve_fit( shifter, r, Phi_data.v)
plt.scatter( r, Phi_data, label='data')

plt.plot( r, Phi_analytic, label='Analytic')
#check out this syntax.  popt is a list, the star * expands the list into arguments.
plt.plot( r, shifter( r, *popt), label='shift %0.1f %0.1f'%(popt[0],popt[1]))
plt.legend(loc=0)
plt.savefig('plots_to_sort/ray.png')
