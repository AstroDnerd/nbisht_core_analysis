
from starter2 import *
import xtra_energy

from scipy.optimize import curve_fit
from scipy import stats
from scipy.ndimage import gaussian_filter
import core_proj_three
reload(core_proj_three)
import other_scrubber
reload(other_scrubber)
#import three_loopers_six as TL
import camera_path

class multipro2():
    def __init__(self,mon):
        self.mon=mon
        self.cores_used=[]
        self.profiles_gas={}
        self.profiles_part={}
        self.ds={}

    def run(self, core_list=None, frame_list=None, timescale=0, get_particles=False, save_sorts=False):
        self.timescale=timescale
        suites_to_use=[0]
        if get_particles:
            suites_to_use=[0,1]
        if core_list is None:
            core_list = self.mon.core_list
        if timescale in [2]:
            Ntimes = 4
        else:
            Ntimes = 6
        for core_id in core_list:
            self.profiles_gas[core_id]={}
            self.profiles_part[core_id]={}
            ms = self.mon.get_ms(core_id, do_central=True, do_ge=True,do_ke=True)
            self.cores_used.append(core_id)

            if self.timescale==0:
                frame_list = self.mon.frames_from_tsung(core_id, [0.0,0.25,0.5,0.75])

                self.titles=[ r'$t=0$', r'$t=0.25 t_{\rm{sing}}$', r'$t=0.5 t_{\rm{sing}}$', r'$t=0.75 t_{\rm{sing}}$',\
                    r'$t=t_{\rm{sing}}$', r'$t=t_{\rm{sung}}$']
            if self.timescale==2:
                frame_list = self.mon.frames_from_tsung(core_id, [0.0, 0.5])
                self.titles=[ r'$t=0$', r'$t=0.5 t_{\rm{sing}}$', r'$t=t_{\rm{sing}}$', r'$t=t_{\rm{sung}}$']
            rm = rainbow_map(len(frame_list))

            for nframe,frame in enumerate(frame_list):
                self.profiles_gas[core_id][frame]={}
                self.profiles_part[core_id][frame]={}
                print("profile on %s c%04d n%04d"%(self.mon.name, core_id, frame))
                ds = self.mon.get_ds(frame)
                
                for suite in suites_to_use:
                    collector={}
                    sph = self.mon.get_sphere(core_id, frame, 'rinf')
                    nf = self.mon.get_frame_index(frame)

                    if suite == 0:
                        dv = sph[YT_cell_volume]
                        RR = sph['radius']
                        DD = sph[YT_density]
                        EG = sph[YT_grav_energy_2]
                        B2 = sph[YT_magnetic_field_strength]

                        ORDER = np.argsort( RR)
                        if 1:
                            if 0:
                                #probably the right thing to do
                                vel = []
                                for axis in 'xyz':
                                    vel.append( sph['velocity_%s'%axis][ORDER][:10].mean())
                            if 1:
                                #the consistent thing to do
                                vel = ds.arr(ms.vcentral[:,nf], 'code_velocity')
                            scrub = other_scrubber.scrubber(sph, reference_velocity = vel)
                            scrub.compute_ke_rel()
                        EK = scrub.ke_rel
                        vt = scrub.vt_rel
                        vr = scrub.vr_rel
                        vr[ np.isnan(vr)]=0
                        vt[ np.isnan(vt)]=0
                    elif suite == 1:
                        dv = ms.cell_volume[:,nf]
                        RR = ms.r[:,nf]
                        DD = ms.density[:,nf]
                        vr = ms.vr_rel[:,nf]
                        vt = ms.vt2_rel[:,nf]**0.5
                        EG = ms.ge[:,nf]
                        EK = ms.ke_rel[:,nf]
                        B2 = this.looper.tr.c([core_id],'magnetic_field_strength')
                        ORDER = np.argsort( RR)

                    rho_sort = DD[ORDER]
                    RR_sort = RR[ORDER]
                    dv_sort = dv[ORDER]
                    b2_sort = B2[ORDER]


                    M_cuml = np.cumsum( DD[ORDER]*dv[ORDER])
                    d2_cuml = np.cumsum( DD[ORDER]**2*dv[ORDER])
                    V_cuml = np.cumsum( dv[ORDER])
                    EG_cuml = np.abs(np.cumsum( EG[ORDER]*dv[ORDER]))/V_cuml
                    EK_cuml = np.cumsum( EK[ORDER]*dv[ORDER])/V_cuml
                    vr_cumsum = np.cumsum( vr[ORDER]*dv_sort)/V_cuml
                    vt_cumsum = np.cumsum( vt[ORDER]*dv_sort)/V_cuml
                    #pdb.set_trace()

                    flux = -( DD[ORDER]*vr[ORDER]*dv[ORDER])
                    r_bins = np.geomspace( 2e-4, 32/128, 32)
                    rcen = 0.5*(r_bins[1:]+r_bins[:-1])
                    digitized = np.digitize( RR_sort, r_bins)
                    mean_flux  =nar([ flux[ digitized == i].sum() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])
                    mass_quant  =nar([ M_cuml[ digitized == i].mean() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])
                    drho  =nar([ rho_sort[ digitized == i].mean() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])

                    #collector['rho']= M_cuml/V_cuml
                    collector['b2_sort']=b2_sort
                    collector['rho']=M_cuml/V_cuml
                    collector['drho']=rho_sort
                    collector['dvr'] = vr[ORDER]
                    collector['dvt'] = vt[ORDER]
                    collector['V_cuml']=V_cuml

                    collector['R']=RR_sort
                    collector['vr_cumsum']=vr_cumsum
                    if np.isnan(vr_cumsum).any():
                        pdb.set_trace()

                    collector['vt_cumsum']=vt_cumsum
                    collector['energy']=EK_cuml/EG_cuml
                    collector['EGmean']=EG_cuml
                    collector['EKmean']=EK_cuml

                    ok = ~np.isnan(mean_flux)
                    collector['rbins']=rcen[ok]
                    collector['mdot']=mean_flux[ok]/mass_quant[ok]*colors.tff
                    collector['mass_quant']=mass_quant[ok]

                    sqrtG=np.sqrt(colors.G)
                    collector['britton']=sqrtG*(collector['rho'])/(np.abs(vr_cumsum.v)+1)

                    if save_sorts:
                        #collector['rho_sort']=colors.density_units*rho_sort
                        collector['rho_sort']=rho_sort
                        collector['dv_sort']=dv_sort
                        collector['vr_sort']=vr[ORDER]
                        collector['vt_sort']=vt[ORDER]

                    if suite==0:
                        self.profiles_gas[core_id][frame]=collector
                    elif suite==1:
                        self.profiles_part[core_id][frame]=collector



