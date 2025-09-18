
from starter2 import *
import convex_hull_tools as CHT

def density(mon,cloud):
    ds = mon.get_ds(0)
    cg = ds.covering_grid(0,[0]*3,[128]*3)


    #min and max, quantized to zone center
    dx = 1/128
    this_min = cloud.min(axis=0)
    this_max = cloud.max(axis=0)
    #this_min = nar([0.5*dx]*3)
    #this_max = nar([1-0.5*dx]*3)
    this_min = (this_min//dx+0.5)*dx
    this_max = (this_max//dx+0.5)*dx

    xyz = np.mgrid[this_min[0]:this_max[0]+dx:dx, 
                      this_min[1]:this_max[1]+dx:dx, 
                      this_min[2]:this_max[2]+dx:dx]
    xyzp = np.column_stack([xyz[0].flatten(),xyz[1].flatten(),xyz[2].flatten()])
    mask_ok = CHT.in_hull(xyzp,cloud)
    yes = xyzp[mask_ok]

    #periodic shift
    yes[ yes < 0 ] += 1
    yes[ yes > 1 ] -= 1
    N = 1/dx
    ind = yes//dx
    index = (ind[:,2]+N*ind[:,1]+N*N*ind[:,0]).astype('int')
    density = cg[YT_density].flatten()[index]
    return density.mean(),density.std()
    #full_x = cg['z'].flatten()
    #print('did I do it right?',(full_x[index].v - yes[:,2]))


    return -1, -1


if 0:
    #ask about points in the 
    ok = ((all_points >= this_min)*(all_points <= this_max)).all(axis=1)
    fast_mask = CHT.in_hull(all_points[ok],cloud)
    mask = copy.copy(ok)
    mask[mask] *= fast_mask

