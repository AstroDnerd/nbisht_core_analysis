from starter2 import *

class scrubber():
    def __init__(self,obj,center=None,reference_velocity=None,do_velocity=True, do_magnetic=False):
        self.obj=obj
        self.scrub(center,reference_velocity,do_velocity=do_velocity, do_magnetic=do_magnetic)
                

    def scrub(self,center=None, reference_velocity=None, axis=0, do_velocity=True, do_magnetic=False):

        self.this_x = self.obj[YT_x]
        self.this_y = self.obj[YT_y]
        self.this_z = self.obj[YT_z]

        self.density = self.obj[YT_density]
        self.cell_volume = self.obj[YT_cell_volume]
        self.mass = self.density*self.cell_volume
        self.mass_total=self.mass.sum(axis=0)

        self.mean_x = np.mean(self.this_x,axis=0)
        self.mean_y = np.mean(self.this_y,axis=0)
        self.mean_z = np.mean(self.this_z,axis=0)
        self.mean_center=nar([self.mean_x,self.mean_y,self.mean_z])
        self.mean_xc= np.sum(self.this_x*self.mass,axis=0)/self.mass_total
        self.mean_yc= np.sum(self.this_y*self.mass,axis=0)/self.mass_total
        self.mean_zc= np.sum(self.this_z*self.mass,axis=0)/self.mass_total
        self.mean_center_density=nar([self.mean_xc,self.mean_yc,self.mean_zc])

        self.rx_rel=self.this_x-self.mean_x
        self.ry_rel=self.this_y-self.mean_y
        self.rz_rel=self.this_z-self.mean_z

        self.rx_relc=self.this_x-self.mean_xc
        self.ry_relc=self.this_y-self.mean_yc
        self.rz_relc=self.this_z-self.mean_zc

        self.r2 = self.rx_rel**2+self.ry_rel**2+self.rz_rel**2
        self.r2c = self.rx_relc**2+self.ry_relc**2+self.rz_relc**2

        self.r=np.sqrt(self.r2)
        self.rc=np.sqrt(self.r2c)

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
            self.mean_vx = self.raw_vx.mean()
            self.mean_vy = self.raw_vy.mean()
            self.mean_vz = self.raw_vz.mean()

            self.mass_mean_vx = self.raw_vx*self.mass/self.mass_total
            self.mass_mean_vy = self.raw_vy*self.mass/self.mass_total
            self.mass_mean_vz = self.raw_vz*self.mass/self.mass_total

            self.rel_vx = self.raw_vx-self.mean_vx
            self.rel_vy = self.raw_vy-self.mean_vy
            self.rel_vz = self.raw_vz-self.mean_vz
            self.rel_vmag = (self.rel_vx**2+self.rel_vy**2+self.rel_vz**2)**(0.5)

            self.rx_hat= self.rx_rel/self.r
            self.ry_hat= self.ry_rel/self.r
            self.rz_hat= self.rz_rel/self.r

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

    def compute_ge(self):
        self.gx = self.obj[YT_grav_x]
        self.gy = self.obj[YT_grav_y]
        self.gz = self.obj[YT_grav_z]
        self.ge = -1/(np.pi*8*colors.G)*(self.gx**2+self.gy**2+self.gz**2)
    def compute_ke(self,core_id):
        self.ke = 0.5*self.density*(self.raw_vx**2+self.raw_vy**2+self.raw_vz**2)
    def compute_ke_rel(self,core_id):
        self.ke_rel = 0.5*self.density*(self.rel_vx**2+self.rel_vy**2+self.rel_vz**2)



