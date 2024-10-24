

from starter2 import *
import track_loader as TL
import colors
#import ucore
#reload(ucore)


import dtools.davetools as dt


sim_list=['u501','u502','u503']
import track_loader as TL
TL.load_tracks(sim_list)
import monster
monster.load(sim_list)

class meanie():
    def __init__(self,mon):
        self.mon=mon
        self.cores_used = []
        self.frames_used = []
        self.means = {}
    def run(self,core_list,sphere_type='rmax',frames='movie'):
        if frames == 'movie':
            frame_slice=slice(None)
        elif frames == 'short':
            frame_slice = np.zeros_like(self.mon.frames,dtype='bool')
            frame_slice[::10]=True
            frame_slice[-1]=True
        elif frames == 'realshort':
            frame_slice=slice(-1,None,10)
            frame_slice = np.zeros_like(self.mon.frames,dtype='bool')
            frame_slice[0]=True;frame_slice[-1]=True

        frame_list=self.mon.frames[frame_slice]
        self.frames_used=frame_list
        self.cores_used = core_list
        self.means['B']={}
        self.means['rho']={}
        for frame in frame_list:
            self.means['B'][frame]=[]
            self.means['rho'][frame]=[]
            for core_id in core_list:
                print('core frame',core_id, frame)
                sph = self.mon.get_sphere(core_id,frame,sphere_type)
                mtot = sph['cell_mass'].sum()
                vtot = sph['cell_volume'].sum()
                self.means['B'][frame].append( (sph['magnetic_field_strength']*sph['cell_mass']).sum()/mtot)
                self.means['rho'][frame].append(mtot/vtot)



def ploot(brho):
    rm = dt.rainbow_map(len(brho.frames_used))
    norm = mpl.colors.Normalize(vmin=0,vmax=1)
    cmap=copy.copy(mpl.cm.get_cmap("jet"))
    cmap.set_over([0.5,0.5,0.5])
    colorbar = mpl.cm.ScalarMappable(norm=norm,cmap=cmap).to_rgba

    fig,ax=plt.subplots(1,1)
    ext_b = dt.extents()
    ext_d = dt.extents()
    for nf,frame in enumerate(brho.frames_used):
        nf = brho.mon.get_frame_index(frame)
        time = brho.mon.all_times[nf]/colors.tff
        t_tsung = [time/mon.get_tsung(core_id) for core_id in brho.cores_used]
        col = [colorbar(ttt) for ttt in t_tsung]
        this_b = brho.means['B'][frame]
        this_d = brho.means['rho'][frame]
        #ax.scatter(this_d,this_b,c=[rm(nf)]*len(this_d))
        ax.scatter(this_d,this_b,c=col)
        ext_b(this_b)
        ext_d(this_d)
    ax.set(xscale='log',yscale='log',xlim=ext_d.minmax,ylim=ext_b.minmax)
    fig.savefig('%s/brho_%s'%(plot_dir,brho.mon.name))
def ploot2(brho,brho2):
    rm = dt.rainbow_map(len(brho.frames_used))

    fig,axes=plt.subplots(1,2)
    ax0=axes[0];ax1=axes[1]

    ext_b = dt.extents()
    ext_d = dt.extents()
    for nf,frame in enumerate(brho.frames_used):
        this_b_0 = brho.means['B'][frame]
        this_d_0 = brho.means['rho'][frame]
        this_b_2 = brho2.means['B'][frame]
        this_d_2 = brho2.means['rho'][frame]
        ax0.scatter(this_d_0,this_d_2,c=[rm(nf)]*len(this_d_0))
        ax1.scatter(this_b_0,this_b_2,c=[rm(nf)]*len(this_d_0))
        ext_b(this_b_0)
        ext_d(this_d_0)
        ext_b(this_b_2)
        ext_d(this_d_2)
    ax0.set(xscale='log',yscale='log',xlim=ext_d.minmax,ylim=ext_d.minmax, xlabel='d r1',ylabel='d rmax')
    ax1.set(xscale='log',yscale='log',xlim=ext_b.minmax,ylim=ext_b.minmax, xlabel='B r1', ylabel='B rmax')
    ax0.plot([1,1e3],[1,1e3],c='k')
    ax1.plot([1,1e3],[1,1e3],c='k')
    fig.tight_layout()
    fig.savefig('%s/brho2_%s'%(plot_dir,brho.mon.name))





if 'brho_r1' not in dir() or True:
    brho_r1={}
    brho_rmax={}
if 'brho_rsmart' not in dir():
    brho_rsmart={}
for sim in sim_list[:1]:
    if sim in brho_rsmart:
        continue
    mon = monster.closet[sim]
    core_list =  mon.this_looper.core_by_mode['A']
    this = meanie(mon)
    this.run(core_list,frames='short',sphere_type='rsmart')
    brho_rsmart[sim]=this
if 0:
    this = meanie(mon)
    this.run(core_list, frames='short',sphere_type='r1')
    brho_r1[sim]=this
    that = meanie(mon)
    that.run(core_list, frames='short',sphere_type='rmax')
    brho_rmax[sim]=that

ploot(brho_rsmart[sim])
#for sim in sim_list[:1]:
#    ploot2(brho_r1[sim], brho_rmax[sim])

