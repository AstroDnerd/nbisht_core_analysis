

from starter2 import *

# PLOT OF THE 2D,3D ALPHAS
fig,ax = plt.subplots(1,1)
a3D = [0.456282,0.479462,0.617871]
a2D_cyl =[0.711603,0.713381,0.731364]
a2D_reg = [0.520261,0.446243,0.496914]
ax.plot(a3D, a2D_cyl, 'bo',alpha=0.8)
ax.plot(a3D, a2D_reg, 'bs',alpha=0.8)
ax.set_xlim(0.2,0.9)
ax.set_ylim(0.2,0.9)
ax.set_xlabel(r'$\alpha_{3D}$')
ax.set_ylabel(r'$\alpha_{2D}$')
fig.savefig('alpha_2D_vs_3D')
plt.close(fig)

