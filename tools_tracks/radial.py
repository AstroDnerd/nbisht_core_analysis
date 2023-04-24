
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

class multipro():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
        self.profiles_gas={}
        self.profiles_part={}

    def run(self, core_list=None, frame_list=None, tsing=None, timescale=0, get_particles=False, save_sorts=False):
        self.timescale=timescale
        this_looper=self.this_looper
        thtr=this_looper.tr
        suites_to_use=[0]
        if get_particles:
            suites_to_use=[0,1]
        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)
        def get_time_index(time):
            index=np.argmin( np.abs( thtr.times/colors.tff-time))
            return index
        Nplots = 4
        if timescale in [2]:
            Ntimes = 4
        else:
            Ntimes = 6
        ext = [extents() for n in range(Nplots+1)]
        for core_id in core_list:
            self.profiles_gas[core_id]={}
            self.profiles_part[core_id]={}
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            ms.get_central_at_once(core_id)
            ms.compute_ge(core_id)
            ms.compute_ke_rel(core_id)
            self.cores_used.append(core_id)

            frame_mask = np.zeros_like(thtr.times, dtype='bool')
            if self.timescale==1:
                self.titles=[]
                for theta in np.linspace(0,1,Ntimes):
                    t = (1-theta)*tsing.tsing_core[core_id]+theta*tsing.tend_core[core_id]
                    self.titles.append(r'%0.2f'%theta)
                    index=get_time_index(t)
                    frame_mask[index]=True


            if self.timescale==2:
                frame_mask[0]=True
                #frame_mask[get_time_index(0.25*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(0.5*tsing.tsing_core[core_id])]=True
                #frame_mask[get_time_index(0.75*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(tsing.tend_core[core_id])]=True
                #half_collapse = 0.5*(tsing.tsing_core[core_id]+tsing.tend_core[core_id])
                #theta=0.5
                #half_collapse = theta*tsing.tsing_core[core_id]+(1-theta)*tsing.tend_core[core_id]
                #frame_mask[get_time_index(half_collapse)]=True
                self.titles=[ r'$t=0$', r'$t=0.5 t_{\rm{sing}}$', r'$t=t_{\rm{sing}}$', r'$t=t_{\rm{sung}}$']
            if self.timescale==0:
                frame_mask[0]=True
                frame_mask[get_time_index(0.25*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(0.5*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(0.75*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(tsing.tend_core[core_id])]=True
                #half_collapse = 0.5*(tsing.tsing_core[core_id]+tsing.tend_core[core_id])
                #theta=0.5
                #half_collapse = theta*tsing.tsing_core[core_id]+(1-theta)*tsing.tend_core[core_id]
                #frame_mask[get_time_index(half_collapse)]=True
                self.titles=[ r'$t=0$', r'$t=0.25 t_{\rm{sing}}$', r'$t=0.5 t_{\rm{sing}}$', r'$t=0.75 t_{\rm{sing}}$',\
                    r'$t=t_{\rm{sing}}$', r'$t=t_{\rm{sung}}$']
            
            frame_list=thtr.frames[frame_mask]
            rm = rainbow_map(len(frame_list))

            for nframe,frame in enumerate(frame_list):
                self.profiles_gas[core_id][frame]={}
                self.profiles_part[core_id][frame]={}
                print("profile on %s c%04d n%04d"%(this_looper.sim_name, core_id, frame))
                ds = this_looper.load(frame)
                xtra_energy.add_energies(ds)
                xtra_energy.add_gdotgradrho(ds)
                nf = np.where( this_looper.tr.frames == frame)[0][0]

                center = nar([ms.mean_xc[nf], ms.mean_yc[nf],ms.mean_zc[nf]])
                msR = ms.rc
                msR[ msR<1/2048]=1/2048
                
                for suite in suites_to_use:
                    collector={}
                    MaxRadius=msR[:,nf].max()
                    Radius = max([8.0/128, MaxRadius])
                    rsph = ds.arr(Radius,'code_length')
                    sph = ds.sphere(center,rsph)

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

                    flux = -( DD[ORDER]*vr[ORDER]*dv[ORDER])
                    r_bins = np.geomspace( 2e-4, 32/128, 32)
                    rcen = 0.5*(r_bins[1:]+r_bins[:-1])
                    digitized = np.digitize( RR_sort, r_bins)
                    mean_flux  =nar([ flux[ digitized == i].sum() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])
                    mass_quant  =nar([ M_cuml[ digitized == i].mean() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])

                    #collector['rho']= M_cuml/V_cuml
                    collector['b2']=b2_sort
                    collector['rho']=colors.density_units* M_cuml/V_cuml
                    collector['V_cuml']=V_cuml

                    collector['R']=RR_sort
                    collector['vr_cumsum']=vr_cumsum

                    collector['vt_cumsum']=vt_cumsum
                    collector['energy']=EK_cuml/EG_cuml
                    collector['EGmean']=EG_cuml
                    collector['EKmean']=EK_cuml

                    ok = ~np.isnan(mean_flux)
                    collector['rbins']=rcen[ok]
                    collector['mdot']=mean_flux[ok]/mass_quant[ok]*colors.tff
                    collector['mass_quant']=mass_quant[ok]

                    if suite==0:
                        self.profiles_gas[core_id][frame]=collector
                    elif suite==1:
                        self.profiles_part[core_id][frame]=collector

                    if save_sorts:
                        #collector['rho_sort']=colors.density_units*rho_sort
                        collector['rho_sort']=rho_sort
                        collector['dv_sort']=dv_sort
                        collector['vr_sort']=vr[ORDER]
                        collector['vt_sort']=vt[ORDER]




