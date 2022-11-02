
from starter2 import *


# PLOT OF THE 2D,3D ALPHAS
fig,ax = plt.subplots(1,1)
a3D = [0.456282,0.479462,0.617871] 
a2D_cyl128 =[0.711603,0.713381,0.731364]
a2D_cyl128frb = [0.62,0.54,0.65]
a2D_cylinf = [0.61,0.55,0.57]
a2D_cylinffrb = [0.63,0.62,0.71]
a2D_can = [a2D_cyl128,a2D_cylinf]
a2D_canfrb = [a2D_cyl128frb,a2D_cylinffrb]
a2D_midcyl128 =[0.59,0.53,0.63]
a2D_midcyl128frb = [0.44,0.34,0.53]
a2D_midcylinf = [0.51,0.44,0.49]
a2D_midcylinffrb = [0.51,0.49, 0.61]
a2D_midcan = [a2D_midcyl128,a2D_midcylinf]
a2D_midcanfrb = [a2D_midcyl128frb,a2D_midcylinffrb] 

#for i in range():
#ax.plot(a2D_cyl128, a3D, 'bo',alpha=0.7)
#ax.plot(a2D_cyl128frb, a3D, 'go',alpha=0.7)
ax.plot(a2D_cylinf, a3D, 'ro',alpha=0.7)
ax.plot(a2D_cylinffrb, a3D, 'mo',alpha=0.7)

#ax.plot(a2D_midcyl128, a3D, 'bs',alpha=0.7)
#ax.plot(a2D_midcyl128frb, a3D, 'gs',alpha=0.7)
ax.plot(a2D_midcylinf, a3D, 'rs',alpha=0.7)
ax.plot(a2D_midcylinffrb, a3D, 'ms',alpha=0.7)
ax.axline((0, 0), slope=1, c='k')

ax.set_xlim(0.2,0.9)
ax.set_ylim(0.2,0.9)
ax.set_xlabel(r'$\alpha_{2D}$')
ax.set_ylabel(r'$\alpha_{3D}$')
fig.savefig('alpha_3D_vs_2D')
plt.close(fig)