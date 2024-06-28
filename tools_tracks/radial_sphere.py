
from starter2 import *
import xtra_energy
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import pcolormesh_helper as pch
import colors
import movie_frames 
G = colors.G
plt.close('all')
import scipy.stats

def smother(arr):
    ok = ~np.isnan(arr)
    ab = arr+0
    ab[~ok]=0
    ab=gaussian_filter(ab,2)
    ab[~ok]=np.nan
    return ab
def dIdt(prof):
    print('run dIdt')

    prof.fields['Ismooth'] = smother( prof.fields['I'])
    prof.fields['dI'] = (prof.fields['Ismooth'][:,2:]-prof.fields['Ismooth'][:,:-2])/prof.dt
    #prof.fields['dI'] = (prof.fields['Ismooth'][:,2:]-prof.fields['Ismooth'][:,:-2])/prof.dt
    #prof.fields['dI'] = smother( prof.fields['dI'])
    prof.fields['dIsmooth'] = smother( prof.fields['dI'])
    #prof.fields['dIsmooth'] = prof.fields['dI']+0
    prof.fields['ddI'] = (prof.fields['dIsmooth'][:,2:]-prof.fields['dIsmooth'][:,:-2])/prof.ddt#/prof.fields['I'][:,2:-2]
    prof.fields['dI'] /= prof.fields['Ismooth'][:,1:-1]/colors.tff
    prof.fields['ddI'] /= prof.fields['Ismooth'][:,2:-2]/colors.tff**2

class profiler():
    def __init__(self,mon,core_id):
        self.mon=mon
        self.core_id=core_id
        self.field_list=['M','particle_count','I','dI','ddI','EK','EG','EB','ET']
        self.fields={}
    def run(self, sphere_type='rmax', frames='movie'):
        if frames == 'movie':
            frame_slice=slice(None)
        elif frames == 'short':
            frame_slice=slice(None,None,10)

        frame_list=self.mon.frames[frame_slice]
        self.times=self.mon.times[frame_slice]#/colors.tff
        self.dt   = 0.5*(self.times[2:]-self.times[:-2])
        self.tcen = 0.5*(self.times[2:]+self.times[:-2])
        self.ddt  = 0.5*(self.tcen[2:]-self.tcen[:-2])
        self.ddtcen  =     (self.times[2:-2])
        core_id=self.core_id
        ms = self.mon.get_ms(core_id, do_central=True, do_ge=True,do_ke=True)
        self.r_bins = np.geomspace( 2e-4, 32/128, 32)
        self.r_cen = 0.5*(self.r_bins[1:]+self.r_bins[:-1])
        for field in self.field_list:
            self.fields[field]=np.zeros([len(self.r_bins)-1,len(frame_list)])
        for nf,frame in enumerate(frame_list):
            print('prof',frame)
            sph = self.mon.get_sphere(core_id,frame,sphere_type)
            scrub = self.mon.scrub_sphere(core_id,frame,sphere_type)
            xyz, r = self.mon.particle_recenter(core_id, frame, sph.center.v)
            dv = scrub.cell_volume
            RR = scrub.r
            DD = scrub.density
            ORDER = np.argsort( RR)
            RR_srt = RR[ORDER]
            DD_srt = DD[ORDER]
            dv_srt = dv[ORDER]

            EG_srt = scrub.ge[ORDER]
            EK_srt = scrub.ke_rel[ORDER]
            EB_srt = sph['magnetic_energy'][ORDER]
            ET_srt = DD_srt*np.log(DD_srt)

            M_srt = DD_srt*dv_srt
            I_srt = DD_srt*RR_srt**2

            M_cuml = np.cumsum( M_srt )
            V_cuml = np.cumsum( dv_srt)
            EG_cuml = np.cumsum( EG_srt*dv_srt)
            EK_cuml = np.cumsum( EK_srt*dv_srt)
            EB_cuml = np.cumsum( EB_srt*dv_srt)
            ET_cuml = np.cumsum( ET_srt*dv_srt)
            I_cuml = np.cumsum(I_srt*dv_srt)

            order_particle = np.argsort(r)
            rp_srt = r[order_particle]
            np_cuml = np.cumsum(np.ones_like(rp_srt))
            particle_count, bins, ind = scipy.stats.binned_statistic(rp_srt, np_cuml, bins=self.r_bins, statistic='mean')
            self.fields['particle_count'][:,nf]=particle_count


            Mbinned, bins, ind = scipy.stats.binned_statistic(RR_srt, M_cuml, bins=self.r_bins,statistic='mean')
            EGbinned, bins, ind = scipy.stats.binned_statistic(RR_srt, EG_cuml, bins=self.r_bins,statistic='mean')
            EKbinned, bins, ind = scipy.stats.binned_statistic(RR_srt, EK_cuml, bins=self.r_bins,statistic='mean')
            EBbinned, bins, ind = scipy.stats.binned_statistic(RR_srt, EB_cuml, bins=self.r_bins,statistic='mean')
            ETbinned, bins, ind = scipy.stats.binned_statistic(RR_srt, ET_cuml, bins=self.r_bins,statistic='mean')
            Ibinned, bins, ind = scipy.stats.binned_statistic(RR_srt, I_cuml, bins=self.r_bins, statistic='mean')
            self.fields['M'][:,nf]=Mbinned
            self.fields['I'][:,nf]=Ibinned
            self.fields['EK'][:,nf]=EKbinned
            self.fields['EG'][:,nf]=EGbinned
            self.fields['EB'][:,nf]=EBbinned
            self.fields['ET'][:,nf]=ETbinned


class surface():
    def __init__(self,mon,core_id):
        self.mon=mon
        self.core_id=core_id
        self.field_list=['SK','SB','ST','I','dI','ddI']
        self.fields={}
    def run(self):
        frame_slice=slice(None)
        #frame_slice=slice(None,None,10)
        frame_list=self.mon.frames[frame_slice]
        self.times=self.mon.times[frame_slice]/colors.tff
        core_id=self.core_id
        ms = self.mon.get_ms(core_id, do_central=True, do_ge=True,do_ke=True)
        self.r_bins = np.geomspace( 2e-4, 32/128, 32)
        self.r_cen = 0.5*(self.r_bins[1:]+self.r_bins[:-1])
        for field in self.field_list:
            self.fields[field]=np.zeros([len(self.r_bins)-1,len(frame_list)])
        for nf,frame in enumerate(frame_list):
            sph = self.mon.get_sphere(core_id,frame,'rinf')
            scrub = self.mon.scrub_sphere(core_id,frame,'rinf')
            dv = scrub.cell_volume
            RR = scrub.r
            DD = scrub.density
            vx = scrub.rel_vx
            vy = scrub.rel_vy
            vz = scrub.rel_vz
            x  = scrub.this_x
            y  = scrub.this_y
            z  = scrub.this_z
            vr = scrub.vr_rel
            ds = scrub.cell_volume**(2/3)

            ORDER = np.argsort( RR)
            RR_srt = RR[ORDER]

            SK = (x*vx+y*vy+z*vz)*vr*ds
            SK_srt = SK[ORDER]

            def binner_chicken_dinner(surf):
                binned, bins, ind = scipy.stats.binned_statistic(RR_srt, surf, bins=self.r_bins,statistic='sum')
                return binned
            self.fields['SK'][:,nf]=binner_chicken_dinner(SK_srt)


