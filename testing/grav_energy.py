from starter2 import *
import xtra_energy
plt.clf()

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
    directory='/data/cb1/Projects/P19_CoreSimulations/new_sims/u25b_little_big_2'
    rsphere=0.25
setname = "%s/DD%04d/data%04d"

frame=1
if 'ds' not in dir():
    ds = yt.load(setname%(directory,frame,frame))

if 0:
    #make some projections
    yt.ProjectionPlot(ds,0,YT_density).save('plots_to_sort/%s_n%04d'%(sim,frame))
    yt.ProjectionPlot(ds,0,YT_potential_field).save('plots_to_sort/%s_n%04d'%(sim,frame))

if 1:

    xtra_energy.add_energies(ds)
    G = ds['GravitationalConstant']/(np.pi*4)

    ad = ds.all_data()
  #------------- this is the del phi squared for my data -------------------------------------------#
    phi_grav_squ = ad['grav_energy']
 #---------------------------------------------------------------------------------------------------#   
    energy_all_space = 0.5*np.sum(phi_grav_squ*ad['cell_volume'])
    print(ad['cell_mass'].sum())
    M = (ad['cell_mass'][ ad['density']>2]).sum()

#-------------------------------- calculating the potential using fourier transform---------------------------------------#
    #rho = ds.covering_grid(1,[0.0,0.0,0.0],[128,128,128],fields = ['density'])
    #rhohat = np.fft.fftn(rho['density'])
    #kx,ky,kz = 2*math.pi*np.mgrid[0:128:1,0:128:1,0:128:1]

    #rx,ry,rz = np.mgrid[0:128:1,0:128:1,0:128:1]/128 -0.5
    #rt = np.sqrt(rx**2+ry**2+rz**2)

    #k_squ = kx**2+ky**2+kz**2
    #phihat = -4*np.pi*G*(rhohat.real/k_squ)
    #phihat[0][0][0] = 0
    #phi = np.fft.ifftn(phihat)
    #phi_real = phi.real
    #plt.imshow(phi_real.sum(axis=1))
    #plt.savefig('plots_to_sort/sum')
if 0:
    #a ray from the origin
    #Use either this or the next one
    foot = (0.5,0.5)
    ray = ds.ortho_ray(0,foot)
    plt.clf()
    xc,yc,zc=[0.5]*3
    #r = np.sqrt((ray[YT_x].v-xc)**2 + (ray[YT_y].v-yc)**2 + (ray[YT_z].v-zc)**2)
    Phi_data=ray['GravPotential']

if 1:
    #all the points in the sphere
    ad = ds.all_data()
    xc,yc,zc=[0.5]*3
    r = np.sqrt( (ad[YT_x].v-xc)**2+(ad[YT_y].v-yc)**2+(ad[YT_z].v-zc)**2)
    Phi_data = ad[YT_potential_field]
    index_data = np.argsort(Phi_data)
    r_sort = r[index_data]
    number_of_bins = 32
    r_bin = np.linspace(np.min(r),np.max(r),number_of_bins)
    array_energy_r = []

#------------------------------------- this is used to calculate the energy inside a radiial bin for data, I used this to find the binding energy per radial bin  ----------------------------#
    for i in range(len(r_bin)):
        mask_i = r < r_bin[i]
        energy_all_space = 0.5*np.sum(phi_grav_squ[mask_i]*ad['cell_volume'][mask_i])
        array_energy_r.append(energy_all_space)

#------------------------- values for analytic potential energy del_phi^2--------#
nabla_phi_squ  = np.zeros_like(r)
ok2 = r>rsphere
angle = -(M**2*G)/(4*np.pi)
nabla_phi_squ[ok2] = angle/(r[ok2]**4)
nabla_phi_squ[~ok2] = angle*(r[~ok2]**2/rsphere**6)
#--------------------- calculate the analytic binding energy per radial bin --------------------------------#
e_analy = np.zeros_like(r_bin)
ok3 = r_bin > rsphere
coe_e_out = (-3*G*M**2)/5
e_analy[ok3] = (0.5*G*M**2)/r_bin[ok3] + coe_e_out/rsphere
coe_e_a = (-0.5*G*M**2)/5
e_analy[~ok3] = coe_e_a*(r_bin[~ok3]**5/rsphere**6)
#---------------------- potential energy data-------------------------------------#
#---------------------------- This section is not needed ------------------------------------------------------------#
#Dan, make sure you understand what I'm doing here.
Phi_analytic = np.zeros_like( r)
ok1 = r>rsphere
Phi_analytic[ok1] = -G*M/r[ ok1]
Phi_analytic[~ok1] = -G*M*(3*rsphere**2-r[~ok1]**2)/(2*rsphere**3)

coeff_energy = -(6*G*M**2)/5
energy_analytic = coeff_energy/rsphere


def shifter(x,b):
    #Shifter slides the analytic solution Phi_analytic around to match the data
    if hasattr(Phi_analytic,'v'):
        Phi = Phi_analytic.v
    else:
        Phi = Phi_analytic
    return Phi + b

from scipy.optimize import curve_fit

#fit the analytic, which is stored in the function shifter
popt,pcov = curve_fit( shifter, r, Phi_data.v)
#----------------------------------------------------------------------------------------------
if 0:
    plt.scatter( r, Phi_data, label='data')

    plt.plot( r, Phi_analytic,c='r', label='Analytic')
#check out this syntax.  popt is a list, the star * expands the list into arguments.
    plt.plot( r, shifter( r, *popt), c='k',label='shift %0.1f'%(popt[0]))
    plt.xlabel('r')
    plt.ylabel(r'$(1/(8 \pi G)\nabla \Phi ^2(r)$')
    plt.legend(loc=0)
    plt.savefig('plots_to_sort/ray_2')

#-------------------- plot nabla phi squ for both data and for the analytic---------------------------------#
if 1:

    plt.scatter( r, phi_grav_squ, label='data')

#plt.scatter(rt.flatten(),phi_real.flatten(), label = 'transform')

    plt.plot( r, nabla_phi_squ, c='r', label='Analytic')
#check out this syntax.  popt is a list, the star * expands the list into arguments.
   # plt.plot( r, shifter( r, *popt), label='shift %0.1f'%(popt[0]))

    plt.xlabel('r')
    plt.ylabel(r'$(8 \pi G)^{-1}(\nabla \Phi(r))^2$')
    plt.legend(loc=0)
    plt.yscale('symlog',linthresh=1)
    plt.xscale('log')
    plt.savefig('plots_to_sort/phi_grav_squ_sphere',bbox_inches='tight')
#--------------------------------------- this is a plot of the binding energy for different radial bin ------------------------------------------------------------#
if 0:

    plt.scatter( r_bin,array_energy_r, label='data')

#plt.scatter(rt.flatten(),phi_real.flatten(), label = 'transform')

    plt.plot( r_bin,e_analy, label='Analytic')
#check out this syntax.  popt is a list, the star * expands the list into arguments.
   # plt.plot( r, shifter( r, *popt), label='shift %0.1f'%(popt[0]))
    plt.legend(loc=0)
    plt.savefig('plots_to_sort/energy_vs_r_sphere')
