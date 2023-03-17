from starter2 import *

class scrubber():
    def __init__(self,obj,reference_velocity=None,do_velocity=True, do_magnetic=False):
        self.obj=obj
        self.scrub(reference_velocity,do_velocity=do_velocity, do_magnetic=do_magnetic)
                

    def scrub(self, reference_velocity=None, axis=0, do_velocity=True, do_magnetic=False):

        self.reference_velocity=reference_velocity

        #this_x has had the center and periodic wrap removed
        self.this_x = self.obj['periodic_x']
        self.this_y = self.obj['periodic_y']
        self.this_z = self.obj['periodic_z']

        self.density = self.obj[YT_density]
        self.cell_volume = self.obj[YT_cell_volume]
        self.mass = self.density*self.cell_volume
        self.mass_total=self.mass.sum(axis=0)



        #alias.  Lazy coding, we've already removed the mean.
        self.rx_rel=self.this_x
        self.ry_rel=self.this_y
        self.rz_rel=self.this_z


        self.r = self.obj['radius']

        self.rx_hat= self.rx_rel/self.r
        self.ry_hat= self.ry_rel/self.r
        self.rz_hat= self.rz_rel/self.r

        if do_magnetic:
            self.bx = self.obj[YT_magnetix_x]
            self.by = self.obj[YT_magnetix_y]
            self.bz = self.obj[YT_magnetix_z]
            self.b2 = self.bx*self.bx+self.by*self.by+self.bz*self.bz

        if do_velocity:
            self.raw_vx = self.obj[YT_velocity_x]
            self.raw_vy = self.obj[YT_velocity_y]
            self.raw_vz = self.obj[YT_velocity_z]

            self.sqr_vx = self.raw_vx**2
            self.sqr_vy = self.raw_vy**2
            self.sqr_vz = self.raw_vz**2

            self.raw_v2 = self.raw_vx**2+self.raw_vy**2+self.raw_vz**2
            #print("KLUDGE: using raw mean for velocity")

            if self.reference_velocity is None:
                self.mean_vx = self.raw_vx.mean()
                self.mean_vy = self.raw_vy.mean()
                self.mean_vz = self.raw_vz.mean()
            else:
                self.mean_vx = self.reference_velocity[0]
                self.mean_vy = self.reference_velocity[1]
                self.mean_vz = self.reference_velocity[2]

            self.mass_mean_vx = self.raw_vx*self.mass/self.mass_total
            self.mass_mean_vy = self.raw_vy*self.mass/self.mass_total
            self.mass_mean_vz = self.raw_vz*self.mass/self.mass_total

            self.rel_vx = self.raw_vx-self.mean_vx
            self.rel_vy = self.raw_vy-self.mean_vy
            self.rel_vz = self.raw_vz-self.mean_vz
            self.rel_vmag = (self.rel_vx**2+self.rel_vy**2+self.rel_vz**2)**(0.5)


            self.norm_r = (self.rx_hat**2+self.ry_hat**2+self.rz_hat**2)**(0.5)

            self.vr_raw = self.rx_hat*self.raw_vx+\
                      self.ry_hat*self.raw_vy+\
                      self.rz_hat*self.raw_vz
            self.vr_rel = self.rx_hat*self.rel_vx+\
                      self.ry_hat*self.rel_vy+\
                      self.rz_hat*self.rel_vz
            self.vr_x = self.vr_rel*self.rx_hat
            self.vr_y = self.vr_rel*self.ry_hat
            self.vr_z = self.vr_rel*self.rz_hat
            self.vt2_raw = (self.raw_vx-self.vr_raw*self.rx_hat)**2+\
                       (self.raw_vy-self.vr_raw*self.ry_hat)**2+\
                       (self.raw_vz-self.vr_raw*self.rz_hat)**2
            self.vt2_rel = (self.rel_vx-self.vr_rel*self.rx_hat)**2+\
                       (self.rel_vy-self.vr_rel*self.ry_hat)**2+\
                       (self.rel_vz-self.vr_rel*self.rz_hat)**2
            self.vt_x = self.rel_vx-self.vr_rel*self.rx_hat
            self.vt_y = self.rel_vy-self.vr_rel*self.ry_hat
            self.vt_z = self.rel_vz-self.vr_rel*self.rz_hat
            self.vt_rel=np.sqrt(self.vt2_rel)

    def compute_ge(self):
        #self.gx = self.obj[YT_grav_x]
        #self.gy = self.obj[YT_grav_y]
        #self.gz = self.obj[YT_grav_z]
        self.gx = self.obj[YT_acceleration_x]
        self.gy = self.obj[YT_acceleration_y]
        self.gz = self.obj[YT_acceleration_z]
        self.ge = -1/(np.pi*8*colors.G)*(self.gx**2+self.gy**2+self.gz**2)
    def compute_g_radial(self):
        #self.gx = self.obj[YT_grav_x]
        #self.gy = self.obj[YT_grav_y]
        #self.gz = self.obj[YT_grav_z]
        self.gx = self.obj[YT_acceleration_x]
        self.gy = self.obj[YT_acceleration_y]
        self.gz = self.obj[YT_acceleration_z]
        self.gr = self.rx_rel*self.gx+self.ry_rel*self.gy+self.rz_rel*self.gz
    def compute_ke(self):
        self.ke = 0.5*self.density*(self.raw_vx**2+self.raw_vy**2+self.raw_vz**2)
    def compute_ke_rel(self):
        self.ke_rel = 0.5*self.density*(self.rel_vx**2+self.rel_vy**2+self.rel_vz**2)



