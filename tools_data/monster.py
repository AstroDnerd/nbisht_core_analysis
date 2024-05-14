
from starter2 import *
import xtra_energy

import track_loader as TL
import trackage

if 'closet' not in dir() or True:
    closet = {}

def load(sims):
    for sim in sims:
        if sim not in closet:
            print('monster',sim)
            m = boo(TL.tracks[sim])
            closet[sim]=m
import tsing
class boo():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.name = this_looper.sim_name
        self.tr = this_looper.tr
        self.frames=self.tr.frames
        self.core_ids = np.unique(self.tr.core_ids)
        self.tsing = tsing.te_tc(this_looper)
        self.tsing.run()
        self.ms = {}
        self.spheres={'r8':{},'r1':{},'rmax':{},'rinf':{}}
        self.r_inflection={}
    def get_r_inflection(self,core_id,frame):
        import r_inflection
        if core_id not in self.r_inflection:
            self.r_inflection[core_id]={}
        if frame not in self.r_inflection[core_id]:
            sph = self.get_sphere(core_id,frame,'r8')
            r_inf = r_inflection.RKEEP(sph)
            self.r_inflection[core_id][frame]=r_inf
        return self.r_inflection[core_id][frame]

    def get_ds(self,frame):
        #don't need to cache, looper does it.
        import xtra_energy
        ds = self.this_looper.load(frame)
        xtra_energy.add_energies(ds)
        xtra_energy.add_gdotgradrho(ds)
        return ds
    def get_ms(self,core_id,**kwargs):
        if core_id in self.ms:
            return self.ms[core_id]
        else:
            this_ms = trackage.mini_scrubber(self.tr,core_id, **kwargs)
            self.ms[core_id]=this_ms
            return this_ms
    def get_frame_index(self,frame):
        nf = np.where( self.tr.frames == frame)[0][0]
        return nf
    def get_time_index(self,time):
        thtr=self.tr
        index=np.argmin( np.abs( thtr.times/colors.tff-time))
        return index
    def get_time_frame(self,time):
        thtr=self.tr
        index=np.argmin( np.abs( thtr.times/colors.tff-time))
        frame = self.frames[index]
        return frame
    def get_tsing(self,core_id):
        return self.tsing.tsing_core[core_id]
    def get_tsung(self,core_id):
        return self.tsing.tend_core[core_id]
    def frames_from_tsung(self,core_id,fracs):
        #returns fractions of tsung and tsing.
        tsung = self.get_tsung(core_id)
        tsing = self.get_tsing(core_id)
        ind=      [self.get_time_frame(tsung*frac) for frac in fracs]
        ind.append(self.get_time_frame(tsing))
        ind.append(self.get_time_frame(tsung))
        ind = np.unique(ind)
        ind.sort()
        return ind

    def get_sphere(self,core_id,frame,sphere_type):
        if core_id not in self.spheres[sphere_type]:
            self.spheres[sphere_type][core_id]={}
        if frame not in self.spheres[sphere_type][core_id]:
            nf = self.get_frame_index(frame)
            ms = self.get_ms(core_id)
            #center = nar([ms.this_x[index,nf], ms.this_y[index,nf],ms.this_z[index,nf]])
            density = ms.density[:,nf]
            if density.max() > 1e2 and True:
                index = np.argmax(density)
                center = nar([ms.this_x[index,nf], ms.this_y[index,nf],ms.this_z[index,nf]])
                geomcenter = nar([ms.mean_xc[nf], ms.mean_yc[nf],ms.mean_zc[nf]])
            else:
                center = nar([ms.mean_xc[nf], ms.mean_yc[nf],ms.mean_zc[nf]])
            #center = nar([ms.mean_xc[nf], ms.mean_yc[nf],ms.mean_zc[nf]])
            ds = self.get_ds(frame)
            if sphere_type == 'r8':
                Radius = 8.0/128
                rsph = ds.arr(Radius,'code_length')
            elif sphere_type == 'r1':
                Radius = 1.0/128
                rsph = ds.arr(Radius,'code_length')
            elif sphere_type == 'rmax':
                msR = ms.rc+0
                msR[ msR<1/2048]=1/2048
                MaxRadius=msR[:,nf].max()
                Radius = max([1.0/128, MaxRadius])
                rsph = ds.arr(Radius,'code_length')
            elif sphere_type == 'rinf':
                rinf = self.get_r_inflection(core_id,frame)
                rsph = ds.arr(rinf,'code_length')

            sph = ds.sphere(center,rsph)
            self.spheres[sphere_type][core_id][frame]=sph
        return self.spheres[sphere_type][core_id][frame]


