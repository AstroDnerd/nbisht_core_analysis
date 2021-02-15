'''
The Deformation Tensor

'''
from starter2 import *
import xtra_energy
import data_locations as dl
reload(dl)
np.set_printoptions(threshold=sys.maxsize)

#import three_loopers as TL   #FIND


# - - - - - - - - - - - - - - -
if 'this_simname' not in dir():
    this_simname = 'u05'

# PICK A CORE
# paper cores of u05: 21,70,85,165,275,297
# TWO SNAPSHOTS, obtained from data_puller.py for core 165
f_0 = 10
f_1 = 11

if 'this_looper' not in dir():
    file_list=glob.glob(dl.every_ten[this_simname])
    this_looper=looper.core_looper(directory=dl.sims[this_simname]) 
    for nfile,fname in enumerate(file_list):
        this_looper.load_loop(fname)
        print( "Reading file %d of %d"%(nfile,len(file_list)))
    thtr = this_looper.tr  
    thtr.sort_time()
    all_cores = np.unique(thtr.core_ids)

for num,i in enumerate(all_cores):
    if i == 165:
        core_id = i

# - - - - - - - - - - - - - - -
'''
if takign xyz is correct as opposed to a more favored field (see bottom), 
we can now choose two particles to represent a pair of points in frame 10 and hence in frame 11
to find da_{i,10}, and dx_{i,11} 
'''
# but wait...n_particles.txt says core 165 has 253 particles, vs len(x_p) tells me there are 779...
# ...perhaps that is because these are actually grid positions.
#x_p = thtr.c([core_id],'x') 
#y_p = thtr.c([core_id],'y')
#z_p = thtr.c([core_id],'z')

# chosen randomly
#a_1 = 600
#a_2 = 200 
#a_3 = 300

'''
the next two attempts throw the error- LinAlgError: Singular matrix
without the [0,0,0], it throws error: LinAlgError: Last 2 dimensions of the array must be square...
Q! how does having a 3 by 3 make sense in terms of an interval...
'''
##da = np.array([ [x_p[a_1][0], y_p[a_1][0]], [x_p[a_2][0], y_p[a_2][0]] ]) 
##dx = np.array([ [x_p[a_1][1], y_p[a_1][1]], [x_p[a_2][1], y_p[a_2][1]] ]) 

##da = np.array([ [x_p[a_1][0], y_p[a_1][0], z_p[a_1][0]], [x_p[a_2][0], y_p[a_2][0], z_p[a_2][0]], [0,0,0] ]) 
##dx = np.array([ [x_p[a_1][1], y_p[a_1][1], z_p[a_1][1]], [x_p[a_2][1], y_p[a_2][1], z_p[a_2][1]], [0,0,0] ]) 

#da = np.array([ [x_p[a_1][0], y_p[a_1][0], z_p[a_1][0]], [x_p[a_2][0], y_p[a_2][0], z_p[a_2][0]], [x_p[a_3][0], y_p[a_3][0], z_p[a_3][0]] ]) 
#dx = np.array([ [x_p[a_1][1], y_p[a_1][1], z_p[a_1][1]], [x_p[a_2][1], y_p[a_2][1], z_p[a_2][1]], [x_p[a_2][1], y_p[a_2][1], z_p[a_2][1]] ]) 

# Code it:
# F = np.linalg.solve(dX[5:],dx[5:])  (dX = dx, dx = da)(that’s confusing)
# I'm under the impression it is the opposite dX = da and dx = dx
#F = np.linalg.solve(da,dx)

# That will make F_i
# Average all those tensors, F'
# of how many tensors.


# - - - - - - - - - - - - - - -
'''
BETTER ATTEMPT:
'''
# chosen randomly from the number of particles in the respecive core
pt_a = 200
pt_b = 400 
pt_c = 600

point_0 = this_looper.snaps[f_0][core_id].pos
point_1 = this_looper.snaps[f_1][core_id].pos

#da = np.array(point_0[pt_a],point_0[pt_b],point_0[pt_c])  #ValueError: only 2 non-keyword arguments accepted
#dx = np.array(point_1[pt_a],point_1[pt_b],point_1[pt_c])

da = [point_0[pt_a],point_0[pt_b],point_0[pt_c]]
dx = [point_1[pt_a],point_1[pt_b],point_1[pt_c]]

F = np.linalg.solve(da,dx)

