from starter2 import *



N = 64*1j
y,x = np.mgrid[0:1:N,0:1:N]

a=1
b=2
vx = 3*(a*x+b*y)
vy = 4*(a*x+b*y)

import xtra_operators as xo
reload(xo)

bx=1
by=2
dds = nar([1./128]*3)
b_dot_grad_vx = bx*xo.gradf2d( vx,0,dds) + by*xo.gradf2d(vx,1,dds)
b_dot_grad_vy = bx*xo.gradf2d( vy,0,dds) + by*xo.gradf2d(vy,1,dds)


fig,ax=plt.subplots(1,1)
ax.set_aspect('equal')
ax.pcolormesh(x,y,vx,shading='nearest')
ax.quiver(0.5,0.5,bx,by)
fig.savefig('plots_to_sort/gradtest.png')
