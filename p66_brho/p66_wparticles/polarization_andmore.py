
'''
# TAKE 2
Fiege et al 2000
''We show in Section 4.2 that helical fields threading a filamentary cloud can result in depolarization
toward the filament axis, otherwise similar in appearance to what would be expected from the above
transverse field scenario. However, the depolarization is due to the 3-dimensional structure of the field in
this case and does not depend on the grain shapes or their alignment.''

polarization_pt_0,1,2.h5
 def _magnetic_field_poloidal(field, data):
     normal = data.get_field_parameter("normal")
     Bfields = ustack(
         [
         data[ftype, "relative_magnetic_field_x"],
         data[ftype, "relative_magnetic_field_y"],
         data[ftype, "relative_magnetic_field_z"],
         ]
     )
     theta = data["index", "spherical_theta"]
     phi = data["index", "spherical_phi"]
     return get_sph_theta_component(Bfields, theta, phi, normal)
 def _magnetic_field_toroidal ...same but w/o theta

# TAKE 1
polarization_0,1,2.h5
polarization into Enzo
'''

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
        self.Q_long = [np.zeros(len(core_list)) for x in range(3)]
        self.Q_short = [np.zeros(len(core_list)) for x in range(3)]
        self.U_long = [np.zeros(len(core_list)) for x in range(3)]
        self.U_short = [np.zeros(len(core_list)) for x in range(3)]

        self.theta_long = [np.zeros(len(core_list)) for x in range(3)]
        self.UQfrac_long = [np.zeros(len(core_list)) for x in range(3)]

        self.theta_short = [np.zeros(len(core_list)) for x in range(3)]
        self.UQfrac_short = [np.zeros(len(core_list)) for x in range(3)]
        self.columnrho_short = [np.zeros(len(core_list)) for x in range(3)]


        self.theta_mean = [np.zeros(len(core_list)) for x in range(3)]
        self.theta_cos = [np.zeros(len(core_list)) for x in range(3)]

        self.columnrho = [np.zeros(len(core_list)) for x in range(3)]
        self.columnrho_cyl = [np.zeros(len(core_list)) for x in range(3)]
        self.columnrho_midcyl = [np.zeros(len(core_list)) for x in range(3)]
        self.dtheta = [np.zeros(len(core_list)) for x in range(3)]
        self.dtantheta = [np.zeros(len(core_list)) for x in range(3)]
        self.dtheta_cyl = [np.zeros(len(core_list)) for x in range(3)]
        self.dtantheta_cyl = [np.zeros(len(core_list)) for x in range(3)]
        self.dtheta_midcyl = [np.zeros(len(core_list)) for x in range(3)]
        self.dtantheta_midcyl = [np.zeros(len(core_list)) for x in range(3)]
        self.dvlos = [np.zeros(len(core_list)) for x in range(3)]
        self.dvlos_cyl = [np.zeros(len(core_list)) for x in range(3)]
        self.dvlos_midcyl = [np.zeros(len(core_list)) for x in range(3)]
        self.rhoave = [np.zeros(len(core_list)) for x in range(3)]
        self.rhoave_cyl = [np.zeros(len(core_list)) for x in range(3)]
        self.rhoave_midcyl = [np.zeros(len(core_list)) for x in range(3)]
        self.bpos = [np.zeros(len(core_list)) for x in range(3)]
        self.bpos_cyl = [np.zeros(len(core_list)) for x in range(3)]
        self.bpos_total = [np.zeros(len(core_list)) for x in range(3)]
        self.bdcf = [np.zeros(len(core_list)) for x in range(3)]
        self.bdcf_cyl = [np.zeros(len(core_list)) for x in range(3)]
        self.bdcf_midcyl = [np.zeros(len(core_list)) for x in range(3)]
        self.bdcf_tancyl = [np.zeros(len(core_list)) for x in range(3)]
        self.bdcf_sm = [np.zeros(len(core_list)) for x in range(3)]
        self.blos = [np.zeros(len(core_list)) for x in range(3)]
        self.blos_cyl = [np.zeros(len(core_list)) for x in range(3)]
        self.blos_midcyl = [np.zeros(len(core_list)) for x in range(3)]
        self.blos_frb = [np.zeros(len(core_list)) for x in range(3)]
        self.btotal = [np.zeros(len(core_list)) for x in range(3)]
        self.btotal_cyl = [np.zeros(len(core_list)) for x in range(3)]
        self.btotal_midcyl = [np.zeros(len(core_list)) for x in range(3)]
        self.btotal_frb = [np.zeros(len(core_list)) for x in range(3)]
        self.bparticles= [np.zeros(len(core_list)) for x in range(3)]


        self.bpoloidal_long = [np.zeros(len(core_list)) for x in range(3)]
        self.bpoloidal_short = [np.zeros(len(core_list)) for x in range(3)]
        self.btoroidal_long = [np.zeros(len(core_list)) for x in range(3)]
        self.btoroidal_short = [np.zeros(len(core_list)) for x in range(3)]

        self.angular_p = [np.zeros(len(core_list)) for x in range(3)]
        self.angular_pmag = np.zeros(len(core_list)) 
        self.angular_phat = [np.zeros(len(core_list)) for x in range(3)]
        self.bpoloidal_sL = np.zeros(len(core_list))
        self.btoroidal_sL = np.zeros(len(core_list)) 
        self.bpoloidal_sL_ave = np.zeros(len(core_list))
        self.btoroidal_sL_ave = np.zeros(len(core_list)) 

        self.bpoloidal_sL_save = np.zeros(len(core_list))
        self.btoroidal_sL_save = np.zeros(len(core_list)) 
        self.btotal_save = np.zeros(len(core_list))
        self.bpolratio = np.zeros(len(core_list))
        self.btorratio = np.zeros(len(core_list))

        self.angular_pmagnitude = np.zeros(len(core_list)) 
        self.angular_pmagnitude_yt = np.zeros(len(core_list)) 


        # THE FRAMES
        the_frame = thtr.frames[-1:] 

        # CORE-LOOP
        for nc,core_id in enumerate(core_list):
            self.cores_used.append(core_id)

            ds = self.this_looper.load(the_frame[0]) 
            #pdb.set_trace() 
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            ms.particle_pos(core_id)

            p49_fields.add_QU(ds)
            the_center = ms.mean_center[:,-1]  #the three coords for the last frame 
            the_normal = [[1,0,0],[0,1,0],[0,0,1]]
            the_radius = 1/128 
            the_area= np.pi * (the_radius**2) 
            
            # MAKE THE OBJECTS:
            xyz = [0,1,2]
            the_cyl = {}
            the_midcyl = {}
            the_sphere = {}

            for i in range(3):
                # if we wanted spheres; no normal
                the_cyl[xyz[i]] = ds.disk(the_center,the_normal[i],the_radius,height=(1,'code_length'))
                the_midcyl[xyz[i]] = ds.disk(the_center,the_normal[i],the_radius,height=the_radius) 
                #pdb.set_trace()
            
            the_sphere = ds.sphere(the_center, the_radius) 
            Q = ['Qx','Qy','Qz']
            U = ['Ux','Uy','Uz']
            L = ['angular_momentum_x','angular_momentum_y','angular_momentum_z']
            B = ['magnetic_field_x','magnetic_field_y','magnetic_field_z']
            # TO OBTAIN mu_o
            #mu = 4*np.pi * 10**-7 H/m --> 1.11x10^-14 s^2/cm^2
            mu = 4*np.pi
            Q_one = 1
            Q_os = 0.5
            Q_liu = 0.28

            # TO OBTAIN BTOTAL - PARTICLES, 3D
            if 1:
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
                            self.bparticles[j][nc] = (b[j] * cell_volume).sum()/cell_volume.sum() 

            for j in range(3): 
                if 1:
                    self.columnrho[j][nc] = (the_sphere['density'] * the_sphere['cell_volume']).sum()/the_area
                    self.columnrho_cyl[j][nc] = (the_cyl[j]['density'] * the_cyl[j]['cell_volume']).sum()/the_area
                    self.columnrho_midcyl[j][nc] = (the_midcyl[j]['density'] * the_midcyl[j]['cell_volume']).sum()/the_area

                    # FRBs TO OBTAIN DISPERSION THETA 
                    # data source: SPHERE
                    pw = yt.ProjectionPlot(ds,j,('gas','density'),center=the_center,width=the_radius,data_source = the_sphere,origin='window')
                    frb = pw.data_source.to_frb(2*the_radius, 512)  #that is...to_frb(width, resolution, center=None)...is resolution okay?
                    # data source: CYLINDER
                    pw_cyl = yt.ProjectionPlot(ds,j,('gas','density'),center=the_center,width=the_radius,data_source = the_cyl[j],origin='window')
                    frb_cyl = pw_cyl.data_source.to_frb(2*the_radius, 512) 
                    # data source: SHORT CYLINDER
                    pw_midcyl = yt.ProjectionPlot(ds,j,('gas','density'),center=the_center,width=the_radius,data_source = the_midcyl[j],origin='window')
                    frb_midcyl = pw_cyl.data_source.to_frb(2*the_radius, 512) 

                    rho = frb['density']
                    rho_cyl = frb_cyl['density']
                    rho_midcyl = frb_midcyl['density']
                    cv = frb['cell_volume']
                    cv_cyl = frb_cyl['cell_volume']
                    cv_midcyl = frb_midcyl['cell_volume']

                    #pdb.set_trace()
                    #plane = [0,1,2]
                    #los = j
                    #plane.pop(los) 
                    #Q_sph0 = frb[Q[plane[0]]] 
                    #U_sph0 = frb[U[plane[0]]]
                    #Q_sph1 = frb[Q[plane[1]]] 
                    #U_sph1 = frb[U[plane[1]]]
                    #Q_sph = Q_sph0 + Q_sph1
                    #U_sph = U_sph0 + U_sph1

                    Q_sph = frb[Q[j]] 
                    U_sph = frb[U[j]]
                    theta = np.arctan2(U_sph,Q_sph) * 180 / np.pi
                    self.dtheta[j][nc] = np.std(theta)   
                    self.dtantheta[j][nc] = np.std(np.tan(theta))   

                    Q_sph_cyl = frb_cyl[Q[j]]  #should just be Q_cyl 
                    U_sph_cyl = frb_cyl[U[j]]
                    theta_cyl = np.arctan2(U_sph_cyl,Q_sph_cyl) * 180 / np.pi
                    self.dtheta_cyl[j][nc] = np.std(theta_cyl)   
                    self.dtantheta_cyl[j][nc] = np.std(np.tan(theta_cyl))   

                    Q_midcyl = frb_midcyl[Q[j]]  
                    U_midcyl = frb_midcyl[U[j]]
                    theta_midcyl = np.arctan2(U_midcyl,Q_midcyl) * 180 / np.pi
                    self.dtheta_midcyl[j][nc] = np.std(theta_midcyl)   
                    self.dtantheta_midcyl[j][nc] = np.std(np.tan(theta_midcyl))   

                    Q_ave = np.mean(Q_sph)
                    U_ave = np.mean(U_sph)
                    self.theta_mean[j][nc] = np.arctan2(U_ave,Q_ave) * 180 / np.pi
                    print('theta mean',self.theta_mean[j][nc])  
                    # COMPARING ANGLES: EDIT MODE
                    if 0:
                        bx = frb['magnetic_field_x']
                        by = frb['magnetic_field_y']
                        bz = frb['magnetic_field_z']
                        bfrb = [bx, by, bz]
                        #pdb.set_trace()  #continue investigating.
                        b_iave = (bfrb[j]*rho*cv).sum()/(rho*cv).sum()
                        b_iave_x = (bx*rho*cv).sum()/(rho*cv).sum()
                        b_iave_y = (by*rho*cv).sum()/(rho*cv).sum()
                        b_iave_z = (bz*rho*cv).sum()/(rho*cv).sum()
                        cos_theta = (bfrb[j] * b_iave)/(bfrb[j] * np.sqrt(b_iave_x**2 + b_iave_y**2 + b_iave_z**2))
                        self.theta_cos[j][nc] = np.arccos(cos_theta) * 180/np.pi
                        print('difference',self.theta_mean[j][nc] - self.theta_cos[j][nc])
 

                    # TO OBTAIN DISPERSION VLOS: frb or ppv edit!!
                    vlos = frb['velocity_%s'%'xyz'[j]] 
                    vave = (vlos*rho).sum()/rho.sum()
                    self.dvlos[j][nc] = np.sqrt(((vlos-vave)**2).sum())  

                    vlos_cyl = frb_cyl['velocity_%s'%'xyz'[j]] 
                    vave_cyl = (vlos_cyl*rho_cyl).sum()/rho_cyl.sum()
                    self.dvlos_cyl[j][nc] = np.sqrt(((vlos_cyl-vave_cyl)**2).sum())  

                    vlos_midcyl = frb_midcyl['velocity_%s'%'xyz'[j]] 
                    vave_midcyl = (vlos_midcyl*rho_midcyl).sum()/rho_midcyl.sum()
                    self.dvlos_midcyl[j][nc] = np.sqrt(((vlos_midcyl-vave_midcyl)**2).sum())  

                    # TO OBTAIN RHO ...Myers 'mean density'
                    self.rhoave[j][nc] = (the_sphere['density']*the_sphere['cell_volume']).sum()/the_sphere['cell_volume'].sum()
                    self.rhoave_cyl[j][nc] = (the_cyl[j]['density']*the_cyl[j]['cell_volume']).sum()/the_cyl[j]['cell_volume'].sum()
                    self.rhoave_midcyl[j][nc] = (the_midcyl[j]['density']*the_midcyl[j]['cell_volume']).sum()/the_midcyl[j]['cell_volume'].sum()


                    # TO OBTAIN BDCF 
                    if 0:  #HAVE NOT IMPOSED A CORRECTION FACTOR YET 
                        self.bdcf[j][nc] = Q_liu * np.sqrt(mu * self.rhoave[j][nc]) * (self.dvlos[j][nc]/self.dtheta[j][nc])  
                    # SMALL ANGLE APPROXIMATION aka only for <25 degrees...(there are only very few)
                    if 0:
                        if self.dtheta[j][nc] < 25:
                            self.bdcf_sm[j][nc] = Q_one * np.sqrt(mu * self.rhoave[j][nc]) * (self.dvlos[j][nc]/self.dtheta[j][nc])   

                    if 1: #may also consider Q_os
                        self.bdcf[j][nc] = Q_one * np.sqrt(mu * self.rhoave[j][nc]) * (self.dvlos[j][nc]/self.dtheta[j][nc])   
                        self.bdcf_cyl[j][nc] = Q_one * np.sqrt(mu * self.rhoave_cyl[j][nc]) * (self.dvlos_cyl[j][nc]/self.dtheta_cyl[j][nc])   
                        self.bdcf_midcyl[j][nc] = Q_one * np.sqrt(mu * self.rhoave_midcyl[j][nc]) * (self.dvlos_midcyl[j][nc]/self.dtheta_midcyl[j][nc])   
                        self.bdcf_tancyl[j][nc] = Q_one * np.sqrt(mu * self.rhoave_cyl[j][nc]) * (self.dvlos_cyl[j][nc]/self.dtantheta_cyl[j][nc])   
                    # TO OBTAIN BPOS 
                    if 1:
                        #self.bpos[j][nc] = (frb[B[j]] * the_area).sum()  #OR
                        plane = ['x','y','z']
                        los = j
                        plane.pop(los) 
                        bhorizontal = frb['magnetic_field_%s'%plane[0]] 
                        bvertical = frb['magnetic_field_%s'%plane[1]]  
                        bhorizontal_cyl = frb_cyl['magnetic_field_%s'%plane[0]] 
                        bvertical_cyl = frb_cyl['magnetic_field_%s'%plane[1]]  

                        plane = ['x','y','z']
                        bthree = frb['magnetic_field_%s'%plane[los]]

                        self.bpos_total[j][nc] = (np.sqrt(bhorizontal**2 + bvertical**2 + bthree**2) * rho *cv).sum()/(rho*cv).sum()  
                        self.bpos[j][nc] = (np.sqrt(bhorizontal**2 + bvertical**2) * rho *cv).sum()/(rho*cv).sum()  
                        self.bpos_cyl[j][nc] = (np.sqrt(bhorizontal_cyl**2 + bvertical_cyl**2) * rho_cyl *cv_cyl).sum()/(rho_cyl*cv_cyl).sum()  
                        #print('bpos_total',self.bpos_total[j][nc])

                    # TO OBTAIN BLOS...in the past we have tried to do frb style, we abandoned it (maybe re-visit to try again?)
                    if 1:
                        self.blos[j][nc] = (the_sphere['density'] * the_sphere[B[j]] * the_sphere['cell_volume']).sum()/the_sphere['gas','cell_mass'].sum()
                        self.blos_cyl[j][nc] = (the_cyl[j]['density'] * the_cyl[j][B[j]] * the_cyl[j]['cell_volume']).sum()/the_cyl[j]['gas','cell_mass'].sum()
                        self.blos_midcyl[j][nc] = (the_midcyl[j]['density'] * the_midcyl[j][B[j]] * the_midcyl[j]['cell_volume']).sum()/the_midcyl[j]['gas','cell_mass'].sum()
                        cm = frb['cell_mass']
                        self.blos_frb[j][nc] = (rho * bthree * cv).sum()/cm.sum()
                    # TO OBTAIN BTOT
                    if 1: 
                        bmagfield = frb['magnetic_field_strength']
                        self.btotal[j][nc] = (the_sphere['density']* the_sphere['magnetic_field_strength'] * the_sphere['cell_volume']).sum()/(the_sphere['density']*the_sphere['cell_volume']).sum() 
                        self.btotal_cyl[j][nc] = (the_cyl[j]['density']* the_cyl[j]['magnetic_field_strength'] * the_cyl[j]['cell_volume']).sum()/(the_cyl[j]['density']*the_cyl[j]['cell_volume']).sum() 
                        self.btotal_midcyl[j][nc] = (the_midcyl[j]['density']* the_midcyl[j]['magnetic_field_strength'] * the_midcyl[j]['cell_volume']).sum()/(the_midcyl[j]['density']*the_midcyl[j]['cell_volume']).sum() 
                        self.btotal_frb[j][nc] = (bmagfield * rho * cv).sum()/(rho*cv).sum()   
                        #print('btot_frb',self.btotal_frb[j][nc])

                if 0:
                    self.Q_long[j][nc] = (the_cyl[j][Q[j]] * the_cyl[j]['cell_volume']).sum()
                    self.Q_short[j][nc] = (the_mid_cyl[j][Q[j]] * the_mid_cyl[j]['cell_volume']).sum()

                    self.U_long[j][nc] = (the_cyl[j][U[j]] * the_cyl[j]['cell_volume']).sum()
                    self.U_short[j][nc] = (the_mid_cyl[j][U[j]] * the_mid_cyl[j]['cell_volume']).sum()

                    Q_long =self.Q_long[j][nc] 
                    U_long =self.U_short[j][nc]
                    Q_short =self.Q_short[j][nc] 
                    U_short =self.U_short[j][nc]

                    self.theta_long[j][nc] = np.arctan2(U_long,Q_long) * 180 / np.pi  
                    self.theta_short[j][nc] = np.arctan2(U_short,Q_short) * 180 / np.pi  

                    self.columnrho[j][nc] = (the_cyl[j]['density'] * the_cyl[j]['cell_volume']).sum()
                    self.columnrho_short[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j]['cell_volume']).sum()

                    N_long = self.columnrho[j][nc]  
                    N_short = self.columnrho_short[j][nc]  

                    self.UQfrac_long[j][nc] = np.sqrt(Q_long**2 + U_long**2)/N_long 
                    self.UQfrac_short[j][nc] = np.sqrt(Q_short**2 + U_short**2)/N_short 
                    self.bpoloidal_long[j][nc] =(the_cyl[j]['magnetic_field_poloidal'] * the_cyl[j]['cell_volume']).sum() 
                    self.btoroidal_long[j][nc] = (the_cyl[j]['magnetic_field_toroidal'] * the_cyl[j]['cell_volume']).sum()
                    self.bpoloidal_short[j][nc] =(the_mid_cyl[j]['magnetic_field_poloidal'] * the_mid_cyl[j]['cell_volume']).sum() 
                    self.btoroidal_short[j][nc] = (the_mid_cyl[j]['magnetic_field_toroidal'] * the_mid_cyl[j]['cell_volume']).sum()

                    # could try instead: self.angular_p[nc][j]
                    self.angular_p[j][nc] = (the_sphere['cell_mass'] * the_sphere[L[j]]).sum() 
                    # note how each component is averaged first then we can get the magnitude
                    self.angular_pmag[nc] = self.angular_pmag[nc] + (((the_sphere['cell_mass'].v * the_sphere[L[j]].v).sum())/the_sphere['cell_mass'].v.sum())**2  


            if 0:
                self.angular_pmagnitude[nc] = np.sqrt(self.angular_pmag[nc]) 
                self.angular_pmagnitude_yt[nc] = (the_sphere['cell_mass']*the_sphere['angular_momentum_magnitude']).sum()/the_sphere['cell_mass'].sum()
                #pdb.set_trace()
                for i in range(3):
                    self.angular_phat[i][nc] = self.angular_p[i][nc]/np.sqrt(self.angular_pmag[nc])

                the_mid_Lcyl = ds.disk(the_center,[self.angular_phat[0][nc], self.angular_phat[1][nc], self.angular_phat[2][nc]],the_radius,height=(1,'code_length'))

                if 0:
                    # poloidal and toroidal in this new frame, which would only work for the short cans, not spheres.. 
                    self.bpoloidal_sL[nc] =(the_mid_Lcyl['magnetic_field_poloidal'] * the_mid_Lcyl['gas','cell_volume']).sum() 
                    self.bpoloidal_sL_ave[nc] =(the_mid_Lcyl['magnetic_field_poloidal'] * the_mid_Lcyl['gas','cell_volume']).sum()/the_mid_Lcyl['gas','cell_volume'].sum() 
                    self.btoroidal_sL[nc] = (the_mid_Lcyl['magnetic_field_toroidal'] * the_mid_Lcyl['gas','cell_volume']).sum()
                    self.btoroidal_sL_ave[nc] = (the_mid_Lcyl['magnetic_field_toroidal'] * the_mid_Lcyl['gas','cell_volume']).sum()/the_mid_Lcyl['gas','cell_volume'].sum()
                # poloidal and toroidal on spheres, meaning not on any particular normal
                self.bpoloidal_sL_save[nc] =(the_sphere['magnetic_field_poloidal'] * the_sphere['gas','cell_volume']).sum()/the_sphere['gas','cell_volume'].sum() 
                self.btoroidal_sL_save[nc] = (the_sphere['magnetic_field_toroidal'] * the_sphere['gas','cell_volume']).sum()/the_sphere['gas','cell_volume'].sum()
                self.btotal_save[nc] = (the_sphere['magnetic_field_strength'] * the_sphere['gas','cell_volume']).sum()/the_sphere['gas','cell_volume'].sum()
                self.bpolratio[nc] = self.bpoloidal_sL_save[nc]/self.btotal_save[nc]  
                self.btorratio[nc] = self.btoroidal_sL_save[nc]/self.btotal_save[nc] 


                L_normal = [self.angular_phat[0][nc], self.angular_phat[1][nc], self.angular_phat[2][nc]]
                # IMAGE DENSITY WITH THIS SAME NORMAL, aka at respective oblique angles 
                # for all cores in 'alone'
                cutting_plane = yt.SlicePlot(ds, L_normal,('gas','density'),center=the_center,width=0.1,data_source=the_sphere)   #64: 0.1, 128: 0.05, 256: 0.025
                cutting_plane.annotate_title("Bpol_frac %.3f, Btol_frac %.3f"%(self.bpolratio[nc],self.btorratio[nc]))
                outname = '64_core_id_%d_%s'%(core_id,sim)
                cutting_plane.save(outname)
                print(outname)

        print("check-in")
        #pdb.set_trace()
        # CACHE TO DISK 
        if 0: 
            hfivename = 'p66_brho/b_losdcftot_%s.h5'%(sim)
            Fptr = h5py.File(hfivename,'w')
            if 0:
                Fptr['theta_long']=self.theta_long
                Fptr['theta_short']=self.theta_short
                Fptr['pol_frac_long']=self.UQfrac_long 
                Fptr['pol_frac_short']=self.UQfrac_short

                Fptr['b_poloidal_long']=self.bpoloidal_long
                Fptr['b_poloidal_short']=self.bpoloidal_short
                Fptr['b_toroidal_long']=self.boroidal_long
                Fptr['b_toroidal_short']=self.btoroidal_short
                
                Fptr['b_poloidal_sL']=self.bpoloidal_sL
                Fptr['b_poloidal_sL_ave']=self.bpoloidal_sL_ave
                Fptr['b_toroidal_sL']=self.btoroidal_sL
                Fptr['b_toroidal_sL_ave']=self.btoroidal_sL_ave

                # I think this will only be useful if I have it for all time frames; also, look at L paper.
                Fptr['L_ave']=self.angular_pmagnitude
                Fptr['L_ave_yt']=self.angular_pmagnitude_yt

            Fptr['core_id']=self.cores_used

            Fptr['B_dcf']=self.bdcf
            Fptr['B_dcf_cyl']=self.bdcf_cyl
            Fptr['B_dcf_midcyl']=self.bdcf_midcyl
            Fptr['B_dcf_tancyl']=self.bdcf_tancyl
            Fptr['B_dcf_sm']=self.bdcf_sm

            Fptr['B_pos']=self.bpos
            Fptr['B_pos_cyl']=self.bpos_cyl
            Fptr['B_postotal']=self.bpos_total

            Fptr['B_los']=self.blos
            Fptr['B_los_cyl']=self.blos_cyl
            Fptr['B_los_midcyl']=self.blos_midcyl
            Fptr['B_losfrb']=self.blos_frb

            Fptr['B_tot']=self.btotal
            Fptr['B_tot_cyl']=self.btotal_cyl
            Fptr['B_tot_midcyl']=self.btotal_midcyl
            Fptr['B_totfrb']=self.btotal_frb
            Fptr['B_particles']=self.bparticles

            Fptr['N']=self.columnrho
            Fptr['N_cyl']=self.columnrho_cyl
            Fptr['N_midcyl']=self.columnrho_midcyl
            #Fptr['theta_uq']=self.theta_mean
            #Fptr['theta_los']=self.theta_cos

            Fptr.close()



# MAIN
import three_loopers_six as TL6
if 'clobber' not in dir():
    clobber=True
if 'scope1' not in dir() or clobber:
    scope1=polarization(TL6.loops['u601'])
    core_list1 = TL6.loops['u601'].core_by_mode['Alone']
if 'scope2' not in dir() or clobber:
    scope2=polarization(TL6.loops['u602'])
    core_list2 = TL6.loops['u602'].core_by_mode['Alone']
if 'scope3' not in dir() or clobber:
    scope3=polarization(TL6.loops['u603'])
    core_list3 = TL6.loops['u603'].core_by_mode['Alone']

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
        hfivename = 'p66_brho/h5files/b_losdcftot_%s.h5'%(nt)  #EDIT
        Fptr = h5py.File(hfivename,'r')
        if 0:      
            theta_long = Fptr['theta_long'][()]
            theta_short = Fptr['theta_short'][()]
            diff_theta = abs(theta_short - theta_long)  # EDIT FOR > 90 BELOW 
            relerr_theta = abs((theta_short - theta_long)/theta_short)  

            UQfrac_long = Fptr['pol_frac_long'][()]
            UQfrac_short = Fptr['pol_frac_short'][()] 
            relerr_UQfrac = abs((UQfrac_short - UQfrac_long)/UQfrac_short) 
            
            blong_poloidal = Fptr['b_poloidal_long'][()]
            blong_toroidal = Fptr['b_toroidal_long'][()]
            blong_pt_ratio = blong_poloidal/blong_toroidal 
            bshort_poloidal = Fptr['b_poloidal_short'][()]
            bshort_toroidal = Fptr['b_toroidal_short'][()]
            bshort_pt_ratio = bshort_poloidal/bshort_toroidal  

            b_pol_sL = Fptr['b_poloidal_sL']
            b_pol_sL_ave = Fptr['b_poloidal_sL_ave']
            b_tor_sL = Fptr['b_toroidal_sL']
            b_tor_sL_ave = Fptr['b_toroidal_sL_ave']
            L_ave = Fptr['L_ave']
            L_ave_yt = Fptr['L_ave_yt']
       
        cores_used = Fptr['core_id']
        
        b_dcf = Fptr['B_dcf']
        b_DCF = np.concatenate((b_dcf[0],b_dcf[1],b_dcf[2]))
        b_dcflog = np.log10(b_DCF)

        b_dcf_cyl = Fptr['B_dcf_cyl']
        b_DCF_cyl = np.concatenate((b_dcf_cyl[0],b_dcf_cyl[1],b_dcf_cyl[2]))
        b_dcfcyllog = np.log10(abs(b_DCF_cyl))
        ok_dcf = b_dcfcyllog > 1

        b_dcf_midcyl = Fptr['B_dcf_midcyl']
        b_DCF_midcyl = np.concatenate((b_dcf_midcyl[0],b_dcf_midcyl[1],b_dcf_midcyl[2]))
        b_dcfmidcyllog = np.log10(abs(b_DCF_midcyl))

        b_dcf_tancyl = Fptr['B_dcf_tancyl']
        b_DCF_tancyl = np.concatenate((b_dcf_tancyl[0],b_dcf_tancyl[1],b_dcf_tancyl[2]))
        b_dcftancyllog = np.log10(abs(b_DCF_tancyl))
        ok = b_dcftancyllog > 1

        b_dcfsm = Fptr['B_dcf_sm']   #mainly 0.0s...

        b_pos = Fptr['B_pos']
        b_POS = np.concatenate((b_pos[0],b_pos[1],b_pos[2]))
        b_poslog = np.log10(b_POS)
        #ok = abs(b_poslog) > 1

        b_pos_cyl = Fptr['B_pos_cyl']
        b_POS_cyl = np.concatenate((b_pos_cyl[0],b_pos_cyl[1],b_pos_cyl[2]))
        b_poscyllog = np.log10(b_POS_cyl)
        #ok = b_poscylog > 1

        b_losfrb = Fptr['B_losfrb']
        b_LOSFRB = np.concatenate((b_losfrb[0],b_losfrb[1],b_losfrb[2]))
        b_losfrblog = np.log10(b_LOSFRB)

        b_los = Fptr['B_los']
        b_LOS = np.concatenate((b_los[0],b_los[1],b_los[2]))
        b_loslog = np.log10(b_LOS)

        b_los_cyl = Fptr['B_los_cyl']
        b_LOS_cyl = np.concatenate((b_los_cyl[0],b_los_cyl[1],b_los_cyl[2]))
        b_loscyllog = np.log10(abs(b_LOS_cyl))
        ok_los = b_loscyllog > 1

        b_los_midcyl = Fptr['B_los_midcyl']
        b_LOS_midcyl = np.concatenate((b_los_midcyl[0],b_los_midcyl[1],b_los_midcyl[2]))
        b_losmidcyllog = np.log10(abs(b_LOS_midcyl))
        ok_midlos = b_losmidcyllog > 1

        b_totfrb = Fptr['B_totfrb']
        b_TOTFRB = np.concatenate((b_totfrb[0],b_totfrb[1],b_totfrb[2]))
        b_totfrblog = np.log10(b_TOTFRB)

        b_tot = Fptr['B_tot']
        b_TOT = np.concatenate((b_tot[0],b_tot[1],b_tot[2]))
        b_totlog = np.log10(abs(b_TOT))
        #ok = btotlog > 1

        b_tot_cyl = Fptr['B_tot_cyl']
        b_TOT_cyl = np.concatenate((b_tot_cyl[0],b_tot_cyl[1],b_tot_cyl[2]))
        b_totcyllog = np.log10(abs(b_TOT))

        b_tot_midcyl = Fptr['B_tot_midcyl']
        b_TOT_midcyl = np.concatenate((b_tot_midcyl[0],b_tot_midcyl[1],b_tot_midcyl[2]))
        b_totmidcyllog = np.log10(abs(b_TOT_midcyl))

        b_par = Fptr['B_particles']

        n_rho = Fptr['N']
        n_RHO = np.concatenate((n_rho[0],n_rho[1],n_rho[2]))
        n_rholog = np.log10(n_RHO)

        n_rho_cyl = Fptr['N_cyl']
        n_RHO_cyl = np.concatenate((n_rho_cyl[0],n_rho_cyl[1],n_rho_cyl[2]))
        n_rhocyllog = np.log10(n_RHO_cyl)

        n_rho_midcyl = Fptr['N_midcyl']
        n_RHO_midcyl = np.concatenate((n_rho_midcyl[0],n_rho_midcyl[1],n_rho_midcyl[2]))
        n_rhomidcyllog = np.log10(n_RHO_midcyl)

        the_bins = 128  #default
        dxyz = ['x','y','z']
        color=['b','g','orange']
        relerr_UQfrac_all = []
        diff_theta_all = []
        
        '''
        index_z = []
        index_dcf = []
        for a in range(3):
            for num_los,B in enumerate(b_los_cyl[a]):
                if b_loscyllog[num_los] <= 1:
                    print('which zeeman core, which direction',cores_used[num_los],a)
                    index_z.append(cores_used[num_los])
            for num_dcf,B in enumerate(b_dcf_cyl[a]):
                if b_dcfcyllog[num_dcf] <= 1:
                    print('which dcf core, which direction',cores_used[num_dcf],a)
                    index_dcf.append(cores_used[num_dcf])
            print('zeeman outliers',index_z)
            print('dcf outliers',index_dcf)
        '''


        for i in range(3):
            index = []
            # CORRECTED FOR THETA DIFF > 90
            if 0:
                for j,theta in enumerate(diff_theta[i]):
                    if theta > 90:
                        #print('theta',theta)
                        minusoneeighty = 180 - theta
                        #print('corrected theta',minusoneeighty)
                        diff_theta[i][j] = minusoneeighty

                # CORRECTED FOR POL FRAC ERROR > 5
                for k,frac in enumerate(relerr_UQfrac[i]):
                    if frac > 1:
                        print('polfrac',frac)
                        print('which core, which direction',cores_used[k],i)
                        index.append(k)

            # HISTOGRAMS
            if 0:  # THETA PDFS & CDFS 
                pdf_long, bins_long = np.histogram(theta_long[i], bins=the_bins, density=True)
                cdf_long = np.cumsum(pdf_long)
                bin_centers_long = 0.5*(bins_long[1:]+bins_long[:-1]) 
                plt.plot(bin_centers_long,pdf_long,c=color[i])

                pdf_short, bins_short = np.histogram(theta_short[i], bins=the_bins, density=True)
                cdf_short = np.cumsum(pdf_short)
                bin_centers_short = 0.5*(bins_short[1:]+bins_short[:-1]) 
                plt.plot(bin_centers_short,pdf_short,c=color[i],linestyle='dashed')            
            if 0:  # POL FRACS LONG & SHORT PDFS & CDFS 
                pdf_long, bins_long = np.histogram(UQfrac_long[i], bins=the_bins, density=True)
                cdf_long = np.cumsum(pdf_long)
                bin_centers_long = 0.5*(bins_long[1:]+bins_long[:-1]) 
                plt.plot(bin_centers_long,cdf_long,c=color[i])#alpha=0.5)

                pdf_short, bins_short = np.histogram(UQfrac_short[i], bins=the_bins, density=True)
                cdf_short = np.cumsum(pdf_short)
                bin_centers_short = 0.5*(bins_short[1:]+bins_short[:-1]) 
                plt.plot(bin_centers_short,cdf_short,c=color[i],linestyle='dashed')#,alpha=0.5)

            # SCATTERS
            if 0:  # THETA 
                plt.scatter(theta_long[i], theta_short[i], alpha=0.7)
                plt.axline((0, 0), slope=1, c='k', linewidth=1.0)
                outname = 'polarization_fraction_relerr_%d_%s'%(i,simnames[nt])
            if 0:  # POLARIZATION FRACS 
                plt.scatter(UQfrac_long[i], UQfrac_short[i], alpha=0.7)
                plt.axline((0, 0), slope=1, c='k', linewidth=1.0)

            # RELATIVE ERRORS, histograms
            if 0:  # THETA 
                pdf, bins = np.histogram(relerr_theta[i], bins=the_bins, density=True)
                cdf= np.cumsum(pdf)
                bin_centers = 0.5*(bins[1:]+bins[:-1]) 
                plt.plot(bin_centers,pdf,c=color[i],alpha=0.7)
            if 0:  # SAVE, LABEL
                plt.xlabel(r'$\theta_%d$ rel error'%i)
                plt.ylabel('PDF')
                outname = 'polarization_theta_relerr_%d_%s'%(i,simnames[nt])
            if 0:  # POLARIZATION FRACS 
                #the_bins = 64 
                pdf, bins = np.histogram(relerr_UQfrac[i], bins=the_bins, density=True)
                cdf= np.cumsum(pdf)
                bin_centers = 0.5*(bins[1:]+bins[:-1]) 
                plt.plot(bin_centers,cdf,c=color[i],alpha=0.7)
            if 0:  # SAVE, LABEL
                plt.xlabel(r'$(Q^2 + U^2)^{1/2}/N$ rel error %d'%i)
                plt.ylabel('PDF')
                outname = 'polarization_fraction_relerr_%d_%s'%(i,simnames[nt])

            # 3MT PLOTS 
            if 0: 
                the_bins = 16 
                #the_weights = np.ones(len(np.delete(diff_theta[i],index)))/len(np.delete(diff_theta[i],index))
                pdf, bins = np.histogram(np.delete(diff_theta[i],index), bins=the_bins, density=True)#, weights=the_weights)
                cdf= np.cumsum(pdf)
                bin_centers = 0.5*(bins[1:]+bins[:-1]) 
                plt.plot(bin_centers,pdf,c=color[i],alpha=0.7)
                #plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
                plt.xlabel(r'$\Delta\theta_{x,y,z}$')
                plt.ylabel('PDF')
                outname ='theta_diff_%d_%s'%(i,simnames[nt])
                #outname ='theta_diff_%s'%(simnames[nt])
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
                plt.scatter(relerr_UQfrac[i],diff_theta[i],c=color[i],alpha=0.7)
                plt.xlabel(r'$(P_{short}-P_{long})/P_{short}$')
                plt.ylabel(r'$\theta_{short:x,y,z} - \theta_{long:x,y,z}$')
                plt.xlim(-0.05,1.05)
                plt.ylim(-2,92)
                outname ='thetadiff_vs_relerrUQfrac_%d_%s'%(i,simnames[nt])
                #outname ='thetadiff_vs_relerrUQfrac_%s'%(simnames[nt])
            if 0:  # OR HEAT MAP 
                the_bins = 32 
                fig, ax = plt.subplots(1,1)
                output = pch.simple_phase(np.delete(relerr_UQfrac[i],index), np.delete(diff_theta[i],index), log=False, ax=ax, nBins=the_bins)
                fig.colorbar(output['plot'])
                ax.set_xlim(-0.05,1.05)
                ax.set_ylim(-2,92)
                outname ='thetadiff_vs_relerrUQfrac_heatmap_%d_%s'%(i,simnames[nt])
                fig.savefig(outname)
                plt.close('all')
            if 0:  # TO PLOT ALL AT ONCE
                relerr_UQfrac_all = np.concatenate((relerr_UQfrac_all,np.delete(relerr_UQfrac[i],index)))
                diff_theta_all = np.concatenate((diff_theta_all,np.delete(diff_theta[i],index)))


            # CURIOUS - be mindful of the outliers
            if 0: 
                the_bins = 16 
                pdfp, binsp = np.histogram(blong_toroidal[i], bins=the_bins, density=True)
                pdfps, binsps = np.histogram(bshort_toroidal[i], bins=the_bins, density=True)
                bin_centersp = 0.5*(binsp[1:]+binsp[:-1]) 
                bin_centersps = 0.5*(binsps[1:]+binsps[:-1]) 
                plt.plot(bin_centersp,pdfp,c=color[i],alpha=0.7)
                plt.plot(bin_centersps,pdfps,c=color[i],alpha=0.7,linestyle='dashed')
                plt.xlabel(r'$B_{\phi}$ LONG SHORT')
                plt.ylabel(r'PDF LONG SHORT')
                #outname ='Bpol_short_PDF_%d_%s'%(i,simnames[nt])
                outname ='Btol_longshort_PDF_%s'%(simnames[nt])
            if 0:
                plt.scatter(bshort_toroidal[i],bshort_poloidal[i],c=color[i],alpha=0.7)
                plt.axline((0, 0), slope=1, c='k', linewidth=0.8)
                plt.xlabel(r'$B_{\phi}$ SHORT')
                plt.ylabel(r'$B_{z}$ SHORT')
                #outname ='Bpol_vs_Btol_short_%d_%s'%(i,simnames[nt])
                outname ='Bpol_vs_Btol_short_%s'%(simnames[nt])
            if 0:
                plt.scatter(bshort_pt_ratio[i],UQfrac_short[i],c=color[i],alpha=0.7)
                plt.xlim(-50,50)
                plt.ylim(-0.05,0.85)
                plt.xlabel(r'$B_{z}/B_{\phi}$ SHORT')
                plt.ylabel(r'$P_{short}$')
                #outname ='Bpoltol_ratio_long_%d_%s'%(i,simnames[nt])
                outname ='Bpoltol_ratioVsP_short_zoom50_%s'%(simnames[nt])
            if 0:
                plt.scatter(bshort_pt_ratio[i],UQfrac_short[i],c=color[i],alpha=0.7)
                plt.xlabel(r'$B_{z}/B_{\phi}$ SHORT')
                plt.ylabel(r'$P_{short}$')
                outname ='Bpoltol_ratio_short_%d_%s'%(i,simnames[nt])
                outname ='Bpoltol_ratio_short_%s'%(i,simnames[nt])



            # SAVE ONE DIRECTION AT A TIME
            if 0:
                plt.savefig(outname)
                print('plotted_%d'%i)
                plt.close('all')


        # ALL DIRECTIONS AT ONCE  

        index_z = []
        for num_los,B in enumerate(b_loscyllog):
            if b_loscyllog[num_los] <= 1:
                index_z.append(num_los)
        #print('zeeman outliers',index_z)

        if 0: #B_LOSPOSTOTS
            #pdb.set_trace()
            if 1:
                pfit = np.polyfit(n_rhomidcyllog[ok_los],b_dcfmidcyllog[ok_los], 1)
                alpha = pfit[0]
                b_poso = pfit[1]
                n_Rho = np.linspace(n_rhomidcyllog[ok_los].min(),n_rhomidcyllog[ok_los].max(),num=len(n_rhomidcyllog[ok_los]))
                N_Rho = 10 ** n_Rho
                B_two = 10 ** (alpha*n_Rho + b_poso)
                plt.plot(N_Rho,B_two,color='k',linestyle='dashed')
            if 1:
                pearX,pearY = scipy.stats.pearsonr(n_RHO_midcyl[ok_los],abs(b_DCF_midcyl)[ok_los])

            #plt.scatter(np.delete(n_RHO_cyl,index_z), np.delete(abs(b_DCF_cyl),index_z),c = 'orange', alpha=0.4)
            plt.scatter(n_RHO_midcyl[ok_los], abs(b_DCF_midcyl)[ok_los],c = 'orange', alpha=0.4)
            #plt.scatter(n_RHO[ok], b_LOS[ok],c = 'g', alpha=0.4)

            #pdb.set_trace()
            #plt.scatter(abs(b_DCF_cyl), abs(b_DCF_midcyl), c = 'orange', alpha=0.4)
            #plt.scatter(abs(b_TOT_midcyl)[ok_los], abs(b_DCF_midcyl)[ok_los], c = 'orange', alpha=0.4)
            #plt.axline((0, 0), slope=1, c='k', linewidth=0.5)

            #plt.scatter(n_rho, b_los,c = 'g', alpha=0.4)
            #plt.scatter(n_rho, b_par,c = 'r', alpha=0.4)
            #plt.scatter(n_rho, b_sph,c = 'k', alpha=0.4)
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel(r'$N$ $(cm^{-2})$')
            #plt.xlabel(r'$B_{cyldcf}$ $(\mu G)$')
            #plt.ylabel(r'$B_{dcf}$ orange, $B_{los}$ green $(\mu G)$')
            plt.ylabel(r'$B_{dcf}$ $(\mu G)$')
            plt.title(r'$\kappa = %f, R=%f$'%(alpha,pearX))
            outname ='b_midcyldcf_zok_%s'%(simnames[nt])

        color=['g','orange']  #b_tor blue, b_pol orange
        bdcf_frac = abs(b_DCF_midcyl)[ok_los]/abs(b_DCF_cyl)[ok_los] 
        bzeem_frac = abs(b_LOS_midcyl)[ok_los]/abs(b_LOS_cyl)[ok_los] 
        if 1: 
            the_bins = 64 #16 
            pdf1, bins1 = np.histogram(bzeem_frac, bins=the_bins, density=True) 
            pdf2, bins2 = np.histogram(bdcf_frac, bins=the_bins, density=True)
            cdf1= np.cumsum(pdf1)
            cdf2= np.cumsum(pdf2)
            bin_centers1 = 0.5*(bins1[1:]+bins1[:-1]) 
            bin_centers2 = 0.5*(bins2[1:]+bins2[:-1]) 
            plt.plot(bin_centers1,pdf1,alpha=0.7,color=color[0])
            plt.plot(bin_centers2,pdf2,alpha=0.7,color=color[1])
            plt.xlim(-0.25,8.5)
            #plt.xlabel(r'$L_{ave},L_{ave_yt}$ SHORT')
            #plt.xlabel(r'$B_{los}/B_{tot}$, $B_{dcf}/B_{tot}$ z:green, dcf:orange')
            plt.xlabel(r'$B_{midcyldcflos}/B_{cyldcflos}$')
            plt.ylabel(r'PDF')
            outname ='B_midcyldcflos_B_cyldcflos_fracs_zokfix64_%s'%(simnames[nt])

        color=['b','orange']  #b_tor blue, b_pol orange
        if 0: 
            color=['b','orange']  #b_tor blue, b_pol orange; MOVE OUTSIDE IF
            the_bins = 16 
            pdf1, bins1 = np.histogram(np.log10(b_tor_sL_ave), bins=the_bins, density=True)  #there are nans if we do the log10s of b_tol, b_pol
            pdf2, bins2 = np.histogram(np.log10(b_pol_sL_ave), bins=the_bins, density=True)
            cdf1= np.cumsum(pdf1)
            cdf2= np.cumsum(pdf2)
            bin_centers1 = 0.5*(bins1[1:]+bins1[:-1]) 
            bin_centers2 = 0.5*(bins2[1:]+bins2[:-1]) 
            plt.plot(bin_centers1,pdf1,alpha=0.7)
            plt.plot(bin_centers2,pdf2,alpha=0.7,color=color[1])
            #plt.xlim(-50,50)
            #plt.xlabel(r'$L_{ave},L_{ave_yt}$ SHORT')
            plt.xlabel(r'$B_{toroidal_Lave},B_{poloidal_Lave}$ SHORT')
            plt.ylabel(r'PDF')
            outname ='Btolandpol_shortL_ave_log10%s'%(simnames[nt])
            #outname ='L_ave_aveyt_log10%s'%(simnames[nt])
        if 0:
            plt.scatter(b_tor_sL,b_pol_sL,alpha=0.7)
            plt.axline((0, 0), slope=1, c='k', linewidth=0.8)
            #plt.xscale('log')
            #plt.yscale('log')
            plt.xlabel(r'$B_{toroidal}$ SL')
            plt.ylabel(r'$B_{poloidal}$ SL')
            #outname ='Bpol_vs_Btol_longlog_%d_%s'%(i,simnames[nt])
            outname ='Bpol_vs_Btol_shortL_%s'%(simnames[nt])
        if 0:
            plt.xlabel(r'$\theta$')
            plt.ylabel('PDF')
            plt.ylim(0,0.40)
            outname = 'polarization_theta_shortlong_xyz_%s'%(simnames[nt])
        if 0:
            plt.xlabel(r'$\theta$ rel error')
            plt.ylabel('PDF')
            plt.xlim(0.0,1.2)
            plt.ylim(-0.2,7.0)
            outname = 'polarization_theta_relerr_zoomed_xyz_%s'%(simnames[nt])

        if 0:
            plt.xlabel(r'$(Q^2 + U^2)^{1/2}/N$')
            plt.ylabel('PDF')
            outname = 'polarization_fraction_CDFshortlong_xyz_%s'%(simnames[nt])
        if 0: 
            plt.xlabel(r'$(Q^2 + U^2)^{1/2}/N$ rel error')
            plt.ylabel('PDF')
            #plt.xlim(0.0,1.2)
            #plt.ylim(-0.2,9.0)
            outname = 'polarization_fraction_relerr_cdf_bins64_xyz_%s'%(simnames[nt])

        if 0: 
            plt.xlabel(r'$\theta$ LONG')
            plt.ylabel(r'$\theta$ SHORT')
            outname = 'pol_theta_short_vs_long_%s'%(simnames[nt])
        if 0: 
            plt.xlabel(r'$(Q^2 + U^2)^{1/2}/N$ LONG')
            plt.ylabel(r'$(Q^2 + U^2)^{1/2}/N$ SHORT')
            plt.ylim(-0.02,0.85)
            plt.xlim(-0.02,0.6)
            outname = 'pol_fraction_short_vs_long_%s'%(simnames[nt])

        # SAVING ALL AT ONCE
        if 1:
            plt.savefig(outname)
            print('plotted_%d'%i)
            plt.close('all')

        if 0:  # HEAT MAP 
            the_bins = 32 
            fig, ax = plt.subplots(1,1)
            output = pch.simple_phase(relerr_UQfrac_all, diff_theta_all, log=False, ax=ax, nBins=the_bins)
            fig.colorbar(output['plot'])
            ax.set_xlim(-0.05,1.05)
            ax.set_ylim(-2,92)
            plt.xlabel(r'$P_{error:x,y,z}$')
            plt.ylabel(r'$\Delta\theta_{x,y,z}$')
            outname ='thetadiff_vs_relerrUQfrac_heatmap_%s'%(simnames[nt])
            fig.savefig(outname)
            plt.close('all')

        # AND CLOSE H5 FILE
        Fptr.close()


