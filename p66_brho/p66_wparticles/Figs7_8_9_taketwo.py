
from starter2 import *
import davetools
reload(davetools)
import p49_fields
reload(p49_fields)
import math
import pcolormesh_helper as pch
from matplotlib.ticker import PercentFormatter
np.set_printoptions(threshold=sys.maxsize)
# --- --- --- --- --- --- ---

class polarization(): 
    def __init__(self,the_loop):
        self.this_looper = the_loop
        self.cores_used = []
        self.bad_cores = []

    def qtyRun(self,sim,core_list=None): 
        print('inside qtyRun')
        thtr = self.this_looper.tr

        # CORES
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        # FIELDS STORED
        self.columnrho_sph = [np.zeros(len(core_list)) for x in range(3)]
        self.columnrho_cyl = [np.zeros(len(core_list)) for x in range(3)]
        self.columnrho_mcyl_one = [np.zeros(len(core_list)) for x in range(3)]
        self.columnrho_mcyl_two = [np.zeros(len(core_list)) for x in range(3)]
        
        self.dtheta_cyl = [np.zeros(len(core_list)) for x in range(3)]
        self.dtheta_mcyl_one = [np.zeros(len(core_list)) for x in range(3)]
        self.dtheta_mcyl_two = [np.zeros(len(core_list)) for x in range(3)]
 
        self.bdcf_cyl = [np.zeros(len(core_list)) for x in range(3)]
        self.bdcf_mcyl_one = [np.zeros(len(core_list)) for x in range(3)]
        self.bdcf_mcyl_two = [np.zeros(len(core_list)) for x in range(3)]
        
        self.bpos_cyl = [np.zeros(len(core_list)) for x in range(3)]
        self.bpos_mcyl_one = [np.zeros(len(core_list)) for x in range(3)]
        self.bpos_mcyl_two = [np.zeros(len(core_list)) for x in range(3)]

        self.blos_sph = [np.zeros(len(core_list)) for x in range(3)]
        self.blos_cyl = [np.zeros(len(core_list)) for x in range(3)]
        self.blos_mcyl_one = [np.zeros(len(core_list)) for x in range(3)]
        self.blos_mcyl_two = [np.zeros(len(core_list)) for x in range(3)]

        self.baxis_cyl = [np.zeros(len(core_list)) for x in range(3)]
        self.baxis_mcyl_one = [np.zeros(len(core_list)) for x in range(3)]
        self.baxis_mcyl_two = [np.zeros(len(core_list)) for x in range(3)]

        self.btotal_sph = [np.zeros(len(core_list)) for x in range(3)]
        self.btotal_cyl = [np.zeros(len(core_list)) for x in range(3)]
        self.btotal_mcyl_one = [np.zeros(len(core_list)) for x in range(3)]
        self.btotal_mcyl_two = [np.zeros(len(core_list)) for x in range(3)]

        self.bparticles= [np.zeros(len(core_list)) for x in range(3)]
        self.rhoave= [np.zeros(len(core_list)) for x in range(3)]
        self.rhoave_sph = [np.zeros(len(core_list)) for x in range(3)]

        # THE FRAMES
        the_frame = thtr.frames[-1:] 

        # CORE-LOOP
        for nc,core_id in enumerate(core_list):
            self.cores_used.append(core_id)

            ds = self.this_looper.load(the_frame[0]) 
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            ms.particle_pos(core_id)

            p49_fields.add_QU(ds)
            the_center = ms.mean_center[:,-1]  #the three coords for the last frame 
            the_normal = [[1,0,0],[0,1,0],[0,0,1]]
            the_radius_one = 1/128  
            the_radius_two = 2/128 
            the_area_one = np.pi * (the_radius_one**2) 
            the_area_two = np.pi * (the_radius_two**2) 
            
            # MAKE THE OBJECTS:
            xyz = [0,1,2]
            the_cyl = {}
            the_sph = {}
            the_midcyl_one = {}
            the_midcyl_two = {}

            for i in range(3):
                the_cyl[xyz[i]] = ds.disk(the_center,the_normal[i],the_radius_one,height=(1,'code_length'))
                the_midcyl_one[xyz[i]] = ds.disk(the_center,the_normal[i],the_radius_one,height=the_radius_one) 
                the_midcyl_two[xyz[i]] = ds.disk(the_center,the_normal[i],the_radius_two,height=the_radius_one) 
            the_sphere = ds.sphere(the_center, the_radius_one) 

            Q = ['Qx','Qy','Qz']
            U = ['Ux','Uy','Uz']
            B = ['magnetic_field_x','magnetic_field_y','magnetic_field_z']
            mu = 4*np.pi
            Q_one = 1
            Q_os = 0.5
            Q_liu = 0.28

            # B PARTICLES, 3D
            if 0:
                all_frames = thtr.frames
                for nf,frame in enumerate(all_frames): 
                    if frame == all_frames[-1]:
                        mask = ms.compute_unique_mask(core_id, dx=1./2048,frame=nf) 
                        cell_volume = thtr.c([core_id],'cell_volume')[mask,nf]
                        density = thtr.c([core_id],'density')[mask,nf]
                        bx = thtr.c([core_id],'magnetic_field_x')[mask,nf]
                        by = thtr.c([core_id],'magnetic_field_y')[mask,nf]
                        bz = thtr.c([core_id],'magnetic_field_z')[mask,nf]                
                        b = [bx, by, bz]
                        for j in range(3):
                            self.bparticles[j][nc] = (density * b[j] * cell_volume).sum()/(density * cell_volume).sum()   #but a heads up that this is for 1 axis! 
                            self.rhoave[j][nc] = (density * cell_volume).sum()/cell_volume.sum()
            for j in range(3): 
                if 1:
                    self.columnrho_sph[j][nc] = (the_sphere['density'] * the_sphere['cell_volume']).sum()/the_area_one
                    self.columnrho_cyl[j][nc] = (the_cyl[j]['density'] * the_cyl[j]['cell_volume']).sum()/the_area_one
                    self.columnrho_mcyl_one[j][nc] = (the_midcyl_one[j]['density'] * the_midcyl_one[j]['cell_volume']).sum()/the_area_one
                    self.columnrho_mcyl_two[j][nc] = (the_midcyl_two[j]['density'] * the_midcyl_two[j]['cell_volume']).sum()/the_area_two
                    
                    '''
                    # FRBS 
                    # data source: CYLINDER
                    pw_cyl = yt.ProjectionPlot(ds,j,('gas','density'),center=the_center,width=the_radius_one,data_source = the_cyl[j],origin='window')
                    frb_cyl = pw_cyl.data_source.to_frb(2*the_radius_one, 512) 
                    # data source: SHORT CYLINDER, RADIUS ONE
                    pw_midcyl_one = yt.ProjectionPlot(ds,j,('gas','density'),center=the_center,width=the_radius_one,data_source = the_midcyl_one[j],origin='window')
                    frb_mcyl_one = pw_midcyl_one.data_source.to_frb(2*the_radius_one, 512) 
                    # data source: SHORT CYLINDER, RADIUS TWO
                    pw_midcyl_two = yt.ProjectionPlot(ds,j,('gas','density'),center=the_center,width=the_radius_two,data_source = the_midcyl_two[j],origin='window')
                    frb_mcyl_two = pw_midcyl_two.data_source.to_frb(2*the_radius_two, 512) 

                    # DENSITY & CELLVOLUME 
                    rho_cyl = frb_cyl['density']
                    rho_mcyl_one = frb_mcyl_one['density']
                    rho_mcyl_two = frb_mcyl_two['density']

                    cv_cyl = frb_cyl['cell_volume']
                    cv_mcyl_one = frb_mcyl_one['cell_volume']
                    cv_mcyl_two = frb_mcyl_two['cell_volume']

                    # DISPERSION THETA
                    Q_cyl = frb_cyl[Q[j]]  
                    U_cyl = frb_cyl[U[j]]
                    theta_cyl = np.arctan2(U_cyl,Q_cyl) * 180 / np.pi
                    self.dtheta_cyl[j][nc] = np.std(theta_cyl)   

                    Q_mcyl_one = frb_mcyl_one[Q[j]]  
                    U_mcyl_one = frb_mcyl_one[U[j]]
                    theta_midcyl_one = np.arctan2(U_mcyl_one,Q_mcyl_one) * 180 / np.pi
                    self.dtheta_mcyl_one[j][nc] = np.std(theta_midcyl_one)   

                    Q_mcyl_two = frb_mcyl_two[Q[j]]  
                    U_mcyl_two = frb_mcyl_two[U[j]]
                    theta_midcyl_two = np.arctan2(U_mcyl_two,Q_mcyl_two) * 180 / np.pi
                    self.dtheta_mcyl_two[j][nc] = np.std(theta_midcyl_two)   

                    # IDEA: dispersion of tan(theta)
                    #self.dtantheta[j][nc] = np.std(np.tan(theta_cyl))   
 
                    # DISPERSION VLOS
                    vlos = frb_cyl['velocity_%s'%'xyz'[j]] 
                    vave = (vlos*rho_cyl).sum()/rho_cyl.sum()
                    dvlos = np.sqrt(((vlos-vave)**2).sum())  

                    vlos_one = frb_mcyl_one['velocity_%s'%'xyz'[j]] 
                    vave_one = (vlos_one*rho_mcyl_one).sum()/rho_mcyl_one.sum()
                    dvlos_one = np.sqrt(((vlos_one-vave_one)**2).sum())  

                    vlos_two = frb_mcyl_two['velocity_%s'%'xyz'[j]] 
                    vave_two = (vlos_two*rho_mcyl_two).sum()/rho_mcyl_two.sum()
                    dvlos_two = np.sqrt(((vlos_two-vave_two)**2).sum())  
                    '''

                    # AVERAGE DENSITY 
                    self.rhoave_sph[j][nc] = (the_sphere['density'] * the_sphere['cell_volume']).sum()/ the_sphere['cell_volume'].sum()
                    '''
                    rhoave_cyl = (the_cyl[j]['density']*the_cyl[j]['cell_volume']).sum()/the_cyl[j]['cell_volume'].sum()
                    rhoave_one = (the_midcyl_one[j]['density']*the_midcyl_one[j]['cell_volume']).sum()/the_midcyl_one[j]['cell_volume'].sum()
                    rhoave_two = (the_midcyl_two[j]['density']*the_midcyl_two[j]['cell_volume']).sum()/the_midcyl_two[j]['cell_volume'].sum()

                    # B DCF 
                    # No correction factor, Q, imposed yet.
                    self.bdcf_cyl[j][nc] = Q_one * np.sqrt(mu * rhoave_cyl) * (dvlos/self.dtheta_cyl[j][nc])   
                    self.bdcf_mcyl_one[j][nc] = Q_one * np.sqrt(mu * rhoave_one) * (dvlos_one/self.dtheta_mcyl_one[j][nc])   
                    self.bdcf_mcyl_two[j][nc] = Q_one * np.sqrt(mu * rhoave_two) * (dvlos_two/self.dtheta_mcyl_two[j][nc])   
                    '''
                    # B POS ; truth (?) to B DCF, compare one plane at a time
                    if 0:
                        plane = ['x','y','z']
                        los = j
                        plane.pop(los) 
                        bhorizontal_cyl = frb_cyl['magnetic_field_%s'%plane[0]] 
                        bvertical_cyl = frb_cyl['magnetic_field_%s'%plane[1]]  
                        bhorizontal_mcyl_one = frb_mcyl_one['magnetic_field_%s'%plane[0]] 
                        bvertical_mcyl_one = frb_mcyl_one['magnetic_field_%s'%plane[1]]  
                        bhorizontal_mcyl_two = frb_mcyl_two['magnetic_field_%s'%plane[0]] 
                        bvertical_mcyl_two = frb_mcyl_two['magnetic_field_%s'%plane[1]]  

                        self.bpos_cyl[j][nc] = (np.sqrt(bhorizontal_cyl**2 + bvertical_cyl**2) * rho_cyl *cv_cyl).sum()/(rho_cyl*cv_cyl).sum()  
                        self.bpos_mcyl_one[j][nc] = (np.sqrt(bhorizontal_mcyl_one**2 + bvertical_mcyl_one**2) * rho_mcyl_one *cv_mcyl_one).sum()/(rho_mcyl_one*cv_mcyl_one).sum()  
                        self.bpos_mcyl_two[j][nc] = (np.sqrt(bhorizontal_mcyl_two**2 + bvertical_mcyl_two**2) * rho_mcyl_two *cv_mcyl_two).sum()/(rho_mcyl_two*cv_mcyl_two).sum()  

                    # B LOS...in the past we have tried to do frb style, we abandoned it.. 
                    self.blos_sph[j][nc] = (the_sphere['density'] * the_sphere[B[j]] * the_sphere['cell_volume']).sum()/the_sphere['gas','cell_mass'].sum()
                    self.blos_cyl[j][nc] = (the_cyl[j]['density'] * the_cyl[j][B[j]] * the_cyl[j]['cell_volume']).sum()/the_cyl[j]['gas','cell_mass'].sum()
                    self.blos_mcyl_one[j][nc] = (the_midcyl_one[j]['density'] * the_midcyl_one[j][B[j]] * the_midcyl_one[j]['cell_volume']).sum()/the_midcyl_one[j]['gas','cell_mass'].sum()
                    self.blos_mcyl_two[j][nc] = (the_midcyl_two[j]['density'] * the_midcyl_two[j][B[j]] * the_midcyl_two[j]['cell_volume']).sum()/the_midcyl_two[j]['gas','cell_mass'].sum()

                    # B AXIS ; truth to B LOS, compare one axis at a time ... but I don't think there is a difference here...
                    if 0: 
                        self.baxis_cyl[j][nc] = (the_cyl[j]['density'] * the_cyl[j][B[j]] * the_cyl[j]['cell_volume']).sum()/the_cyl[j]['gas','cell_mass'].sum()
                        self.baxis_mcyl_one[j][nc] = (the_midcyl_one[j]['density'] * the_midcyl_one[j][B[j]] * the_midcyl_one[j]['cell_volume']).sum()/the_midcyl_one[j]['gas','cell_mass'].sum()
                        self.baxis_mcyl_two[j][nc] = (the_midcyl_two[j]['density'] * the_midcyl_two[j][B[j]] * the_midcyl_two[j]['cell_volume']).sum()/the_midcyl_two[j]['gas','cell_mass'].sum()

                    # B TOT
                    self.btotal_sph[j][nc] = (the_sphere['density']* the_sphere['magnetic_field_strength'] * the_sphere['cell_volume']).sum()/the_sphere['cell_mass'].sum() 
                    self.btotal_cyl[j][nc] = (the_cyl[j]['density']* the_cyl[j]['magnetic_field_strength'] * the_cyl[j]['cell_volume']).sum()/(the_cyl[j]['density']*the_cyl[j]['cell_volume']).sum() 
                    self.btotal_mcyl_one[j][nc] = (the_midcyl_one[j]['density']* the_midcyl_one[j]['magnetic_field_strength'] * the_midcyl_one[j]['cell_volume']).sum()/(the_midcyl_one[j]['density']*the_midcyl_one[j]['cell_volume']).sum() 
                    self.btotal_mcyl_two[j][nc] = (the_midcyl_two[j]['density']* the_midcyl_two[j]['magnetic_field_strength'] * the_midcyl_two[j]['cell_volume']).sum()/(the_midcyl_two[j]['density']*the_midcyl_two[j]['cell_volume']).sum() 

        print("check-in")
        # CACHE TO DISK 
        if 0: 
            hfivename = 'p66_brho/b_sphere_%s.h5'%(sim)
            Fptr = h5py.File(hfivename,'w')

            #Fptr['B_parts']=self.bparticles
            #Fptr['Rho_avg_parts']=self.rhoave
            Fptr['Rho_avg_sph']=self.rhoave_sph

            Fptr['core_id']=self.cores_used
            '''
            Fptr['Dtheta_cyl']=self.dtheta_cyl
            Fptr['Dtheta_mcyl_one']=self.dtheta_mcyl_one
            Fptr['Dtheta_mcyl_two']=self.dtheta_mcyl_two 

            Fptr['B_dcf_cyl']=self.bdcf_cyl
            Fptr['B_dcf_mcyl_one']=self.bdcf_mcyl_one
            Fptr['B_dcf_mcyl_two']=self.bdcf_mcyl_two

            Fptr['B_pos_cyl']=self.bpos_cyl
            Fptr['B_pos_mcyl_one']=self.bpos_mcyl_one
            Fptr['B_pos_mcyl_two']=self.bpos_mcyl_two
            '''
    
            Fptr['B_los_sph']=self.blos_sph
            #Fptr['B_los_cyl']=self.blos_cyl
            #Fptr['B_los_mcyl_one']=self.blos_mcyl_one
            #Fptr['B_los_mcyl_two']=self.blos_mcyl_two

            #Fptr['B_axis_cyl']=self.baxis_cyl
            #Fptr['B_axis_mcyl_one']=self.baxis_mcyl_one
            #Fptr['B_axis_mcyl_two']=self.baxis_mcyl_two

            Fptr['B_tot_sph']=self.btotal_sph
            '''
            Fptr['B_tot_cyl']=self.btotal_cyl
            Fptr['B_tot_mcyl_one']=self.btotal_mcyl_one 
            Fptr['B_tot_mcyl_two']=self.btotal_mcyl_two 
            '''

            Fptr['N_sph']=self.columnrho_sph
            '''
            Fptr['N_cyl']=self.columnrho_cyl
            Fptr['N_mcyl_one']=self.columnrho_mcyl_one
            Fptr['N_mcyl_two']=self.columnrho_mcyl_two
            '''

            Fptr.close()


# MAIN
import three_loopers_six as TL6
if 'clobber' not in dir():
    clobber=True
if 'scope1' not in dir() or clobber:
    scope1=polarization(TL6.loops['u601'])
    #core_list1 = TL6.loops['u601'].core_by_mode['Alone']
if 'scope2' not in dir() or clobber:
    scope2=polarization(TL6.loops['u602'])
    #core_list2 = TL6.loops['u602'].core_by_mode['Alone']
if 'scope3' not in dir() or clobber:
    scope3=polarization(TL6.loops['u603'])
    #core_list3 = TL6.loops['u603'].core_by_mode['Alone']

simnames = ['u601','u602', 'u603']


for nt,tool in enumerate([scope1,scope2,scope3]):

    # WHICH CORES & RUN IF NOT CACHED
    if 0:
        all_cores = np.unique(tool.this_looper.tr.core_ids)
        core_list = all_cores
        #core_list = all_cores[:1]  #DEBUG
        #core_list = core_list1[:3] 
        tool.qtyRun(nt,core_list=core_list) 

    # ONCE CACHED, HISTOGRAMS SCATTERS & REL ERROR REGARDING POLARIZATION 
    if 1: 
        hfivename = 'p66_brho/h5files/b_sphere_%s.h5'%(nt)  #EDIT
        Fptr = h5py.File(hfivename,'r')

        cores_used = Fptr['core_id']

        b_los_sph = Fptr['B_los_sph']
        b_LOS_sph = np.concatenate((b_los_sph[0],b_los_sph[1],b_los_sph[2]))
        b_lossphlog = np.log10(abs(b_LOS_sph))

        n_sph = Fptr['N_sph']
        N_sph = np.concatenate((n_sph[0],n_sph[1],n_sph[2]))
        n_sphlog = np.log10(N_sph)

        b_tot_sph = Fptr['B_tot_sph']
        b_TOT_sph = np.concatenate((b_tot_sph[0],b_tot_sph[1],b_tot_sph[2]))
        b_totsphlog = np.log10(abs(b_TOT_sph))
        
        rho_sph = Fptr['Rho_avg_sph']
        Rho_sph = np.concatenate((rho_sph[0],rho_sph[1],rho_sph[2]))
        rho_sphlog = np.log10(Rho_sph)
        
        '''
        b_dcf_cyl = Fptr['B_dcf_cyl']
        b_DCF_cyl = np.concatenate((b_dcf_cyl[0],b_dcf_cyl[1],b_dcf_cyl[2]))
        b_dcfcyllog = np.log10(abs(b_DCF_cyl))
        ok_dcf_cyl = b_dcfcyllog > 1

        b_dcf_mcyl_one = Fptr['B_dcf_mcyl_one']
        b_DCF_mcyl_one = np.concatenate((b_dcf_mcyl_one[0],b_dcf_mcyl_one[1],b_dcf_mcyl_one[2]))
        b_dcfmcyllog_one = np.log10(abs(b_DCF_mcyl_one))

        b_pos_mcyl_one = Fptr['B_pos_mcyl_one']
        b_POS_mcyl_one = np.concatenate((b_pos_mcyl_one[0],b_pos_mcyl_one[1],b_pos_mcyl_one[2]))
        b_posmcyllog_one = np.log10(b_POS_mcyl_one)

        b_los_cyl = Fptr['B_los_cyl']
        b_LOS_cyl = np.concatenate((b_los_cyl[0],b_los_cyl[1],b_los_cyl[2]))
        b_loscyllog = np.log10(abs(b_LOS_cyl))
        ok_los = b_loscyllog > 1

        b_los_mcyl_one = Fptr['B_los_mcyl_one']
        b_LOS_mcyl_one = np.concatenate((b_los_mcyl_one[0],b_los_mcyl_one[1],b_los_mcyl_one[2]))
        b_losmcyllog_one = np.log10(abs(b_LOS_mcyl_one))

        b_los_mcyl_two = Fptr['B_los_mcyl_two']
        b_LOS_mcyl_two = np.concatenate((b_los_mcyl_two[0],b_los_mcyl_two[1],b_los_mcyl_two[2]))
        b_losmcyllog_two = np.log10(abs(b_LOS_mcyl_two))

        b_tot_cyl = Fptr['B_tot_cyl']
        b_TOT_cyl = np.concatenate((b_tot_cyl[0],b_tot_cyl[1],b_tot_cyl[2]))
        b_totcyllog = np.log10(abs(b_TOT_cyl))

        b_tot_mcyl_one = Fptr['B_tot_mcyl_one']
        b_TOT_mcyl_one = np.concatenate((b_tot_mcyl_one[0],b_tot_mcyl_one[1],b_tot_mcyl_one[2]))
        b_totmcyllog_one = np.log10(abs(b_TOT_mcyl_one))

        b_tot_mcyl_two = Fptr['B_tot_mcyl_two']
        b_TOT_mcyl_two = np.concatenate((b_tot_mcyl_two[0],b_tot_mcyl_two[1],b_tot_mcyl_two[2]))
        b_totmcyllog_two = np.log10(abs(b_TOT_mcyl_two))

        n_rho_cyl = Fptr['N_cyl']
        n_RHO_cyl = np.concatenate((n_rho_cyl[0],n_rho_cyl[1],n_rho_cyl[2]))
        n_rhocyllog = np.log10(n_RHO_cyl)

        n_rho_mcyl_one = Fptr['N_mcyl_one']
        n_RHO_mcyl_one = np.concatenate((n_rho_mcyl_one[0],n_rho_mcyl_one[1],n_rho_mcyl_one[2]))
        n_rhomcyllog_one = np.log10(n_RHO_mcyl_one)
        n_rho_mcyl_two = Fptr['N_mcyl_two']
        n_RHO_mcyl_two = np.concatenate((n_rho_mcyl_two[0],n_rho_mcyl_two[1],n_rho_mcyl_two[2]))
        n_rhomcyllog_two = np.log10(n_RHO_mcyl_two)

        dtheta_cyl = Fptr['Dtheta_cyl']
        dTHETA_cyl = np.concatenate((dtheta_cyl[0],dtheta_cyl[1],dtheta_cyl[2]))

        dtheta_mcyl_one = Fptr['Dtheta_mcyl_one']
        dTHETA_mcyl_one = np.concatenate((dtheta_mcyl_one[0],dtheta_mcyl_one[1],dtheta_mcyl_one[2]))

        dtheta_mcyl_two = Fptr['Dtheta_mcyl_two']
        dTHETA_mcyl_two = np.concatenate((dtheta_mcyl_two[0],dtheta_mcyl_two[1],dtheta_mcyl_two[2]))
        
        b_par = Fptr['B_parts']
        b_PAR = np.concatenate((b_par[0],b_par[1],b_par[2]))
        b_parlog = np.log10(abs(b_PAR))

        rho_avg = Fptr['Rho_avg_parts']
        rho_AVE = np.concatenate((rho_avg[0],rho_avg[1],rho_avg[2]))
        rho_avelog = np.log10(abs(rho_AVE))
        '''


        the_bins = 128  #default
        color=['b','g','orange']
        # ONE DIRECTION AT A TIME 
        for i in range(3):
            # 3MT PLOTS 
            if 0: 
                the_bins = 16 
                #the_weights = np.ones(len(np.delete(relerr_UQfrac[i],index)))/len(np.delete(relerr_UQfrac[i],index))
                pdf, bins = np.histogram(np.delete(relerr_UQfrac[i],index), bins=the_bins, density=True)#, weights =the_weights)
                cdf= np.cumsum(pdf)
                bin_centers = 0.5*(bins[1:]+bins[:-1]) 
                plt.plot(bin_centers,pdf,c=color[i],alpha=0.7)
                #plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
                plt.xlabel(r'$P_{error: x, y, z}$')
                plt.ylabel('PDF')
                outname = 'UQfrac_relerr_PDF_%d_%s'%(i,simnames[nt])
                #outname = 'UQfrac_relerr_%s'%(simnames[nt])

            if 0: 
                plt.scatter(n_rho_mcyl_one[i],abs(b_dcf_mcyl_one[i]),c='orange',alpha=0.4)
                plt.scatter(n_rho_mcyl_one[i],abs(b_pos_mcyl_one[i]),c='k',alpha=0.4)
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel(r'$N$ $(cm^{-2})$')
                plt.ylabel(r'$B_{dcf}$, $B_{pos}$ $(\mu G)$')
                outname ='B_dcfpos_vsN_%d_%s'%(i,simnames[nt])

            # SAVE ONE DIRECTION AT A TIME
            if 0:
                plt.savefig(outname)
                print('plotted_%d'%i)
                plt.close('all')


        # ALL DIRECTIONS AT ONCE  
        if 1: 
            if 1:
                pfit = np.polyfit(rho_sphlog,b_totsphlog, 1)
                alpha = pfit[0]
                b_poso = pfit[1]
                n_Rho = np.linspace(rho_sphlog.min(),rho_sphlog.max(),num=len(rho_sphlog))
                N_Rho = 10 ** n_Rho
                B_two = 10 ** (alpha*n_Rho + b_poso)
                plt.plot(N_Rho,B_two,color='k',linestyle='dashed')
            if 1:
                pearX,pearY = scipy.stats.pearsonr(Rho_sph,abs(b_TOT_sph))
            if 0:
                pfit2 = np.polyfit(n_sphlog,b_lossphlog, 1)
                alpha2 = pfit2[0]
                b_poso2 = pfit2[1]
                n_Rho2 = np.linspace(n_sphlog.min(),n_sphlog.max(),num=len(n_sphlog))
                N_Rho2 = 10 ** n_Rho2
                B_two2 = 10 ** (alpha2*n_Rho2 + b_poso2)
                plt.plot(N_Rho2,B_two2,color='g',linestyle='dashed')
            if 0:
                pearX2,pearY2 = scipy.stats.pearsonr(N_sph,abs(b_LOS_sph))

            #plt.scatter(n_RHO_mcyl_one, abs(b_LOS_mcyl_one),c = 'k', alpha=0.4)
            plt.scatter(Rho_sph, abs(b_TOT_sph),c = 'k', alpha=0.4)
            #plt.scatter(rho_AVE, abs(b_PAR),c = 'k', alpha=0.4)

            #plt.scatter(abs(b_DCF_cyl), abs(b_DCF_midcyl), c = 'orange', alpha=0.4)
            #plt.scatter(abs(b_TOT_midcyl)[ok_los], abs(b_DCF_midcyl)[ok_los], c = 'orange', alpha=0.4)
            #plt.axline((0, 0), slope=1, c='k', linewidth=0.5)

            plt.xscale('log')
            plt.yscale('log')
            #plt.xlabel(r'$N_{sph,los}$ $(cm^{-2})$')
            plt.xlabel(r'$rho_{ave}$ $(cm^{-3})$')
            plt.ylabel(r'$B_{sph,tot}$ $(\mu G)$')
            #plt.ylabel(r'$B_{particles}$ $(\mu G)$')
            #plt.title(r'$\kappa_{1,2} = %f,%f$ $R_{1,2}=%f,%f$'%(alpha,alpha2,pearX,pearX2))
            plt.title(r'$\kappa= %f$ $R = %f$'%(alpha,pearX))
            outname ='b_sphtot_nook_%s'%(simnames[nt])

        #color=['g','orange']  #b_tor blue, b_pol orange
        #blosm1_frac = abs(b_LOS_mcyl_one)/abs(b_TOT_mcyl_one) 
        #blosm2_frac = abs(b_LOS_cyl)/abs(b_TOT_cyl)
        if 0: 
            the_bins = 64 #16 
            pdf1, bins1 = np.histogram(blosm1_frac, bins=the_bins, density=True) 
            pdf2, bins2 = np.histogram(blosm2_frac, bins=the_bins, density=True)
            #pdf3, bins3 = np.histogram(dTHETA_cyl, bins=the_bins, density=True)
            cdf1= np.cumsum(pdf1)
            cdf2= np.cumsum(pdf2)
            #cdf3= np.cumsum(pdf3)
            bin_centers1 = 0.5*(bins1[1:]+bins1[:-1]) 
            bin_centers2 = 0.5*(bins2[1:]+bins2[:-1]) 
            #bin_centers3 = 0.5*(bins3[1:]+bins3[:-1]) 
            plt.plot(bin_centers1,pdf1,alpha=0.7,color='k')
            plt.plot(bin_centers2,pdf2,alpha=0.7,color='g')
            #plt.plot(bin_centers3,pdf3,alpha=0.7,color='k')
            plt.xlabel(r'$B_{los}/B{tot}$ mcyl1(black), cyl(green)')
            plt.ylabel(r'PDF')
            outname ='B_losVStotfrac_mcyl1_cyl_%s'%(simnames[nt])

        # SAVING ALL AT ONCE
        if 1:
            plt.savefig(outname)
            print('plotted_%d'%i)
            plt.close('all')

        # AND CLOSE H5 FILE
        Fptr.close()


