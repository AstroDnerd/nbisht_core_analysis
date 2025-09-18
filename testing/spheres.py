
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
if 0:
    sim = 'u25'; 
    directory='/scratch1/dcollins/Paper19/u25_little_big'
    directory='/data/cb1/Projects/P19_CoreSimulations/new_sims/u25_little_big'
    rsphere=0.25
if 1:
    sim = 'u25'; 
    directory='/scratch1/dcollins/Paper19/u25_little_big'
    directory='/data/cb1/Projects/P19_CoreSimulations/new_sims/u25b_little_big_2'
    rsphere=0.25
setname = "%s/DD%04d/data%04d"

frame=1
ds = yt.load(setname%(directory,frame,frame))
import xtra_energy
xtra_energy.add_energies(ds)

if 0:
    #make some projections
    yt.ProjectionPlot(ds,0,YT_density).save('plots_to_sort/%s_n%04d'%(sim,frame))
    yt.ProjectionPlot(ds,0,YT_gravpotential).save('plots_to_sort/%s_n%04d'%(sim,frame))

    ad = ds.all_data()
    print(ad['cell_mass'].sum())

    M = (ad['cell_mass'][ ad['density']>2]).sum()

if 1:
    #all the points in the sphere
    ad = ds.all_data()
    xc,yc,zc=[0.5]*3
    r = np.sqrt( (ad[YT_x].v-xc)**2+(ad[YT_y].v-yc)**2+(ad[YT_z].v-zc)**2)
    Phi_data = ad[YT_potential_field]


#analytic result
Phi_analytic = np.zeros_like( r)
ok1 = r>rsphere
Phi_analytic[ok1] = -G*M/r[ ok1]
Phi_analytic[~ok1] = -G*M*(3*rsphere**2-r[~ok1]**2)/(2*rsphere**3)

print('1')
def shifter(x,m,b):
    #Shifter slides the analytic solution Phi_analytic around to match the data
    if hasattr(Phi_analytic,'v'):
        Phi = Phi_analytic.v
    else:
        Phi = Phi_analytic
    return Phi*m+b

if 0:
    #plot the potential with a shift.
    from scipy.optimize import curve_fit

#fit the analytic, which is stored in the function shifter
    popt,pcov = curve_fit( shifter, r, Phi_data.v)
    plt.scatter( r, Phi_data, label='data')

    plt.plot( r, Phi_analytic, label='Analytic')
#check out this syntax.  popt is a list, the star * expands the list into arguments.
    plt.plot( r, shifter( r, *popt), label='shift %0.1f %0.1f'%(popt[0],popt[1]))
    plt.legend(loc=0)
    plt.savefig('plots_to_sort/ray.png')

print('2')
grad_phi_2 = np.zeros_like(r)
grad_phi_2[ok1] = (1/r[ok1]**2)**2
grad_phi_2[~ok1]= (1*r[~ok1]/rsphere**3)**2
grad_phi_2 *= -G*M*M/(4*np.pi)

nabla_phi_squ  = np.zeros_like(r)
ok2 = r>rsphere
angle = -(M**2*G)/(4*np.pi)
nabla_phi_squ[ok2] = angle/(r[ok2]**4)
nabla_phi_squ[~ok2] = angle*(r[~ok2]**2/rsphere**6)
print('3')

grad_phi_2_data = ad[YT_grav_energy]
#mult = grad_phi_2_data.max()/grad_phi_2.max()
fig,ax=plt.subplots(1,1)
sl=slice(None,None,10)
ax.scatter( r[sl], grad_phi_2_data[sl],c='k')
ax.plot( r[sl], grad_phi_2[sl],c='r')
ax.set_yscale('symlog',linthresh=1)
print('5')
fig.savefig('plots_to_sort/grad_phi_2.png')
print('4')
