
from starter2 import *

dx=0.1
x,y,z = np.mgrid[0:1:dx,0:1:dx,0:1:dx]
rr = np.sqrt(x**2+y**2+z**2)

rho = np.zeros_like(rr)
rho[rr>0] = rr[rr>0]**(-2)


the_x =np.log10(rr[rr>0].flatten())
the_y =np.log10(rho[rr>0].flatten())
r,p=scipy.stats.pearsonr(the_x,the_y)
print(r,p)

plt.clf()
plt.scatter(the_x,the_y)
plt.savefig('plots_to_sort/pearson_test.png')
