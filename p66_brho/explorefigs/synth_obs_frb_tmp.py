
'''
synthetic observations, version 4
'''

from cmath import log
from re import X
from turtle import width
from starter2 import *
import annotate_particles_4_cpy
reload(annotate_particles_4_cpy)


class telescope(): 
    def __init__(self,the_loop):
        self.this_looper = the_loop
        self.cores_used = []

    def ZoneFrb(self,sim,core_list=None): 
        print('inside ZoneFrb')
        thtr = self.this_looper.tr

        # CORE_LIST
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        self.synthRho = [np.zeros(len(core_list)) for x in range(3)]
        self.synthField = [np.zeros(len(core_list)) for x in range(3)]
        
        self.synthRho_frb = [np.zeros(len(core_list)) for x in range(3)]
        self.synthField_frb = [np.zeros(len(core_list)) for x in range(3)]

        # THE FINAL FRAME 
        the_frame = thtr.frames[-1:] 

        # CORES
        for nc,core_id in enumerate(core_list):
            ds = self.this_looper.load(the_frame[0]) 
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            ms.particle_pos(core_id)
            self.ms = ms

            # PIECES FOR THE OBJECTS
            the_center = ms.mean_center[:,-1]  #the three coords for the last frame 
            the_normal = [[1,0,0],[0,1,0],[0,0,1]]

            twoeight = 1/128
            the_radius = twoeight
            the_area= np.pi * (the_radius**2) 

            # MAKE THE OBJECTS:
            xyz = [0,1,2]
            xyz = [0]
            the_cyl = {}
            for i in xyz:
                the_cyl[i] = ds.disk(the_center,the_normal[i],the_radius,height=(1,'code_length'))

            # THE DENSITY & FIELD, ZONE METHOD:  2D 
            B = ['magnetic_field_x','magnetic_field_y','magnetic_field_z']
            for j in xyz:
                #self.synthRho[j][nc] = (the_cyl[j]['density'] * the_cyl[j]['cell_volume']).sum()/the_area
                self.synthRho[j][nc] = (the_cyl[j]['density'] * the_cyl[j]['cell_volume']).sum()
                M1 = (the_cyl[j]['density'] * the_cyl[j]['cell_volume']).sum()/ the_cyl[j]['cell_volume'].sum()
                cc = the_cyl[j]
                self.synthField[j][nc] = (the_cyl[j]['density'] * the_cyl[j][B[j]] * the_cyl[j]['cell_volume']).sum()/the_cyl[j]['gas','cell_mass'].sum()

            print('one core density',self.synthRho[0][nc])
            print('one core field',self.synthField[0][nc])

            # THE DENSITY & FIELD, FRB: 2D
            # DENSITY 
            projs_cyl_rho = []
            frbs_cyl_rho = []
            for k in xyz:
                projs_cyl_rho.append(yt.ProjectionPlot(ds,k,('gas','density'),center=the_center,\
                                        width=the_radius, data_source = the_cyl[k], weight_field='cell_volume')) 
                frbs_cyl_rho.append(projs_cyl_rho[k].frb)
                
                length = len(frbs_cyl_rho[k]['gas','density'])  
                self.synthRho_frb[k][nc] = (frbs_cyl_rho[k]['gas','density']).sum()/length
                M2 = (frbs_cyl_rho[k]['gas','density']).sum()
                thefrb = frbs_cyl_rho[k]
                
            # FIELD 
            projs_cyl_B = []
            frbs_cyl_B = []
            for m in xyz:
                projs_cyl_B.append(yt.ProjectionPlot(ds,m,B[m],weight_field =('gas','density'),center=the_center,\
                                        width=the_radius, data_source = the_cyl[m]))
                frbs_cyl_B.append(projs_cyl_B[m].frb)
                
                length = len(frbs_cyl_B[m][B[m]])  
                self.synthField_frb[m][nc] = (frbs_cyl_B[m][B[m]]).sum()

            print('one core density frb',self.synthRho_frb[0][nc])
            print('one core field frb',self.synthField_frb[0][nc])
            print(M1)
            print(M2)
            print('%0.2e'%(M1/(M2/800**2)))
            pdb.set_trace()

# MAIN
import three_loopers_six as TL6
if 'clobber' not in dir():
    clobber=True
if 'scope1' not in dir() or clobber:
    scope1=telescope(TL6.loops['u601'])
if 'scope2' not in dir() or clobber:
    scope2=telescope(TL6.loops['u602'])
if 'scope3' not in dir() or clobber:
    scope3=telescope(TL6.loops['u603'])

simnames = ['u601','u602', 'u603']
for nt,tool in enumerate([scope1]):#,scope2,scope3]): 
    # WHICH CORES
    all_cores = np.unique(tool.this_looper.tr.core_ids)
    core_list = all_cores[2:3]  #DEBUG
    #core_list = all_cores

    # RUN
    tool.ZoneFrb(nt,core_list=core_list) 


