from starter2 import *


dx = 1/1024
x = np.geomspace(dx,1+dx,int(1//dx))
a = 1.86
D = 1e5
rho0 = 1
fig,ax=plt.subplots(1,1)
rho = rho0*(1-x**2)**(-a)
rhodot = 2*a*x*(1-x**2)**(-a-1)
rhodotdot = 2*a*rho0*(1-x**2)**(-a-2)*((2*a+1)*x**2+1)
#ax.plot(x, (1-x**2)**(a+1)/x)
ax.plot(x,rho,c='k')
ax.plot(x,rhodot,c='b')
ax.plot(x,rhodotdot,c='g')
ax.set( yscale='log')
fig.savefig('plots_to_sort/rhodot')
