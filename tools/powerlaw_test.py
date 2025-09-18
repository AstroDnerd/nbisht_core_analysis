from starter2 import *

dx=1/16
x,y,z=np.mgrid[-1:1:dx,-1:1:dx,-1:1:dx]
R = np.sqrt(x**2+y**2+z**2).flatten()
R = R[R>0]
A = 102
B = -2.3
if 1:
    #powerlaw
    Q = (A*R**(2*B+2)).flatten() #+ 10**np.random.random(R.size)
if 0:
    #top hat
    Q = R*0+11
    Rsphere=0.5
    Q[R>Rsphere]=1e-5


plt.close('all')

fig,ax=plt.subplots(1,1)

ax.scatter( R, Q)
axbonk(ax,xscale='log',yscale='log',xlabel='R',ylabel='Q')
ok = R.flatten()>0

def powerlaw(x, alpha, norm):
    return (2*alpha+2)*x+norm
R0=1
popt,pcov=curve_fit( powerlaw, np.log10(R[ok]), np.log10(Q[ok]))
A = 10**popt[1]
alpha=popt[0]
r_line = np.linspace( rmin,rmax,128)
dr = r_line[1:]-r_line[:-1]
rr = 0.5*(r_line[1:]+r_line[:-1])
fit_line=10**( (2*alpha+2)*np.log10(rr)+popt[1])
ax.plot( rr, fit_line,c='g')

#ax.plot( R[ok], 10**powerlaw( np.log10(R[ok]), *popt),c='r')
fig.savefig('plots_to_sort/powerlaw_test.png')
mass1 = (Q*dx**3).sum()
mass2= 4*np.pi*10**popt[1]/(2*popt[0]+5)

mass3 = 4*np.pi*A/( (R0**(2*alpha+4)*(2*alpha+5)))*(rmax**(2*alpha+5)-rmin**(2*alpha+5))
mass5 = (4*np.pi*rr**2*fit_line*dr).sum()

mass4 = 4*np.pi/3*Rsphere**3*Q.max()


print(mass1,mass2,mass3,mass4,mass5)
