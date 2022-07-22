from starter2 import *




def gaussian(the_x,norm,x0,sigma):
    #return norm*np.exp( -(the_x-x0)**2/(2*sigma**2))
    return norm/np.sqrt(2*np.pi*sigma**2)*np.exp( -(the_x-x0)**2/(2*sigma**2))

rho_edge = np.geomspace(1e-10,1e10, 2048)
rho = 0.5*(rho_edge[1:]+rho_edge[:-1])
drho = rho_edge[1:]-rho_edge[:-1]

mu,sig=-0.95,2.3
ggg = gaussian( np.log(rho), 1,mu,sig)
plt.clf()
plt.plot(rho, ggg)
plt.xscale('log')
plt.savefig('plots_to_sort/derp2.png')
A = (ggg*drho/rho).sum()
print('A',A)
mu =  (np.log(rho)*ggg*drho/rho).sum()
print('mu',mu)
var = (( (np.log(rho)-mu)**2*ggg*drho/rho).sum())**0.5
print('sigma',var)


if 0:
    s_edge=np.log(rho_edge)
    s_edge = np.linspace(-100,100,2048)
    s = 0.5*(s_edge[1:]+s_edge[:-1])
    ds = s_edge[1:]-s_edge[:-1]

    a = gaussian( s, 1,2.5 ,3)
    plt.clf()
    plt.plot( s, a)
    plt.savefig('plots_to_sort/derp.png')

    N=(a*ds).sum()
    mu=( (s*a*ds).sum()/N)
    sigma=np.sqrt(( (s-mu)**2*a*ds).sum()/N)
    print(sigma)


