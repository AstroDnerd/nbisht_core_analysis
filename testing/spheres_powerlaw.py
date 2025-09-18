from starter2 import *
import xtra_energy
frame=1
setname = "/data/cb1/Projects/P19_CoreSimulations/new_sims/GravityTestSphere/DD%04d/data%04d"%(frame,frame)
#setname = "/data/cb1/Projects/P19_CoreSimulations/new_sims/u25b_little_big_2/DD%04d/data%04d"%(frame,frame)

ds = yt.load(setname)
G = ds['GravitationalConstant']/(4*np.pi)

xtra_energy.add_energies(ds)

ad = ds.sphere([0.49,0.5,0.5],1)

if 0:
    proj=ds.proj('density',0)
    proj.to_pw().save('plots_to_sort/powerlaw')
    proj=ds.proj(YT_grav_energy,0)
    proj.to_pw().save('plots_to_sort/powerlaw')

R0 = 0.1
if 1:
    R = ad['radius']
    D = ad['density']
    GE= ad[YT_grav_energy]

    RR = R[R<R0]
    DD = D[R<R0]
    GG = np.abs(GE[R<R0])

    rbins = np.linspace(RR.min(),RR.max(),64)

    fitd = np.polyfit( np.log10(RR/R0), np.log10(DD),1)
    fitG = np.polyfit( np.log10(RR/R0), np.log10(GG),1)


fig,ax=plt.subplots(1,1)

ax.scatter( ad['radius'], ad[YT_grav_energy])
phi_line_fit=-10**( fitG[0]*np.log10(RR/R0)+fitG[1])
ax.plot( RR, phi_line_fit,c='g')
ax.set_yscale('symlog',linthresh=1)
ax.set_xscale('log')


axx = ax.twinx()
axx.scatter(ad['radius'], ad['density'],c='r')
#dline=10**( fitd[0]*np.log10(RR/R0)+fitd[1])
rho0=10**(fitd[1])
alpha=fitd[0]
dline = rho0*(RR/R0)**alpha
axx.plot( RR, dline,c='k')
ok = ad[YT_density]>2
M = (ad[YT_cell_mass][ok]).sum()
phi1=-1/(4*np.pi)*G*M*M*RR**(2*alpha+2)/R0**(2*alpha+6)
ax.plot(RR,phi1,c='b')
dv = ad[YT_cell_volume][ok]
ge = ad[YT_grav_energy][ok]
rmin=RR.min()
#E1 = (-G*M*M/R0*(2*alpha+5)).v+(G*M*M/R0**(2*alpha+6)*rmin**(2*alpha+5)).v
E1 = (-G*M*M/R0**(2*alpha+6)*(2*alpha+5))*(R0**(2*alpha+5)-rmin.v**(2*alpha+5))
E2 = (dv*ge).sum()
print(E1/E2)


if 0:
    alpha=fitd[0]
    M = (ad['cell_mass'][ ad['density']>2]).sum()
    rsphere = 0.25
    nabla_phi_squ  = np.zeros_like(rbins)
    ok2 = rbins>rsphere
    angle = -(M**2*G)/(4*np.pi)
    expo1=2*alpha+2
    nabla_phi_squ[ok2] = angle/(rbins[ok2]**(expo1))
    expo2=2*alpha+6
    nabla_phi_squ[~ok2] = angle*(rbins[~ok2]**2/rsphere**(expo2))
    ax.plot( rbins, nabla_phi_squ, c='c')


axbonk(axx,xscale='log', yscale='log')
fig.savefig('plots_to_sort/gradphi.png')
