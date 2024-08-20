

from starter2 import *
import track_loader as TL
import colors
#import ucore
#reload(ucore)
import dtools.vis.pcolormesh_helper as pch

import dtools.davetools as dt


sim_list=['u501','u502','u503']
#sim_list=['u501','u502','u503']
sim_list=['u501','u503']
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
        self.vars={}
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
        #frame_slice[0]=False

        frame_list=self.mon.frames[frame_slice]
        self.frames_used=frame_list
        self.cores_used = core_list
        self.means['B']={}
        self.means['rho']={}
        self.means['rhorho']={}
        self.vars['B']={}
        self.vars['rho']={}
        for frame in frame_list:
            self.means['B'][frame]=[]
            self.means['rho'][frame]=[]
            self.means['rhorho'][frame]=[]
            self.vars['B'][frame]=[]
            self.vars['rho'][frame]=[]
            for core_id in core_list:
                print('core frame',core_id, frame)
                sph = self.mon.get_sphere(core_id,frame,sphere_type)
                mtot = sph['cell_mass'].sum()
                vtot = sph['cell_volume'].sum()
                meanB = (sph['magnetic_field_strength']*sph['cell_mass']).sum()/mtot
                self.means['B'][frame].append(meanB)
                self.vars['B'][frame].append( ((sph['magnetic_field_strength']-meanB)**2*sph['cell_mass']).sum()/mtot)
                meanD=mtot/vtot
                self.means['rho'][frame].append(meanD)
                self.vars['rho'][frame].append( ((sph['density']-meanD)**2*sph['cell_volume'])/vtot)
                self.means['rhorho'][frame].append( (sph['density']*sph['cell_mass']).sum()/mtot)



def ploot(brho,fname):
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
        this_d = brho.means['rhorho'][frame]
        #ax.scatter(this_d,this_b,c=[rm(nf)]*len(this_d))
        ax.scatter(this_d,this_b,c=col)
        ext_b(this_b)
        ext_d(this_d)
    ax.set(xscale='log',yscale='log',xlim=ext_d.minmax,ylim=ext_b.minmax,xlabel='rho',ylabel='B')
    fig.savefig('%s/%s'%(plot_dir,fname))

def ploot2(brho,brho2,vert,horz):
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
    ax0.set(xscale='log',yscale='log',xlim=ext_d.minmax,ylim=ext_d.minmax, xlabel='d r%s'%horz,ylabel='d r%s'%vert)
    ax1.set(xscale='log',yscale='log',xlim=ext_b.minmax,ylim=ext_b.minmax, xlabel='B r%s'%horz, ylabel='B r%s'%vert)
    ax0.plot([1,1e3],[1,1e3],c='k')
    ax1.plot([1,1e3],[1,1e3],c='k')
    fig.tight_layout()
    fig.savefig('%s/bd_vs_bd_%s_%s_%s'%(plot_dir,brho.mon.name,vert,horz))

if 'crap' not in dir():
    crap = {}
def ploot_bckgrnd(brho,fname):
    rm = dt.rainbow_map(len(brho.frames_used))
    norm = mpl.colors.Normalize(vmin=0,vmax=1)
    cmap=copy.copy(mpl.cm.get_cmap("jet"))
    cmap.set_over([0.5,0.5,0.5])
    colorbar = mpl.cm.ScalarMappable(norm=norm,cmap=cmap).to_rgba

    fig,ax=plt.subplots(1,1)
    ext_b = dt.extents()
    ext_d = dt.extents()
    nframes = len(brho.frames_used)
    ncores = len(brho.cores_used)
    bbb = np.zeros([nframes,ncores])
    rrr = np.zeros([nframes,ncores])
    
    tsing = [mon.get_tsing(core_id) for core_id in brho.cores_used]
    times=mon.times/colors.tff+0
    times.shape=times.size,1
    tsing=np.array(tsing)
    #times-tsing
    tsing_index = np.argmin(np.abs((times-tsing)),axis=0)
    #pdb.set_trace()

    print('wtf1')
    target=mon.this_looper.target_frame
    ad = mon.get_ds(target).all_data()
    ad0 = mon.get_ds(0).all_data()


    print('wtf2')
    if 'ball' not in crap:
        crap['ball']=ad['magnetic_field_strength']
    if 'dall' not in crap:
        crap['dall']=ad['density']
    if 'ball0' not in crap:
        crap['ball0']=ad0['magnetic_field_strength']
    if 'dall0' not in crap:
        crap['dall0']=ad0['density']
    ball=crap['ball']
    dall=crap['dall']
    ball0=crap['ball0']
    dall0=crap['dall0']

    print('wtf3')
    #pch.simple_phase(dall,ball,[64,64],ax)
    xbins = np.geomspace(dall.min(),dall.max(),65)
    ybins = np.geomspace(ball.min(),ball.max(),65)
    h2d_accum = np.zeros([64,64])
    print('wtf4')
    for core_id  in brho.cores_used:
        print('burp,',core_id)
        sph = mon.get_sphere(core_id, target, 'rsmart_2')
        balln = sph['magnetic_field_strength']
        dalln = sph['density']
        h2d, binsx, binsy = np.histogram2d(dalln,balln,bins=[xbins,ybins])
        h2d_accum += h2d
    xx,yy= np.meshgrid(xbins,ybins,indexing='ij')
    #ax.pcolormesh(xx,yy,h2d_accum)
    pch.helper(h2d_accum,xbins,ybins,ax=ax)

    h2d0, binsx, binsy = np.histogram2d(dall0,ball0,bins=[xbins,ybins])
    pch.contour(h2d0,xbins,ybins,ax=ax,levels=[1e-3])
    h2dT, binsx, binsy = np.histogram2d(dall,ball,bins=[xbins,ybins])
    pch.contour(h2dT,xbins,ybins,ax=ax,levels=[1e-3])

    #ad = mon.get_ds(target).all_data()
    #ball = ad['magnetic_field_strength']
    #dall = ad['density']
    #pch.simple_phase(dall,ball,[64,64],ax)


    def fixme(arr):
        a = np.array([arr[i] for i in range(len(arr))])
        return a

    for nframe,frame in enumerate(brho.frames_used):
        this_b =fixme(brho.means['B'][frame][:])
        this_d =fixme(brho.means['rhorho'][frame][:])
        bbb[nframe,:] = this_b
        rrr[nframe,:] = this_d
    ext_b(bbb)
    ext_d(rrr)
    for nc in range(ncores):
        core_id = brho.cores_used[nc]
        t_tsung = [time/mon.get_tsung(core_id) for time in mon.times/colors.tff]
        col = [colorbar(ttt) for ttt in t_tsung]
        this_b =bbb[:,nc]#/bbb[tsing_index[nc],nc]
        this_d =rrr[:,nc]#/rrr[tsing_index[nc],nc]
        #ax.scatter(this_d,this_b,c=[rm(nf)]*len(this_d))
        ax.scatter(this_d,this_b,c=col[:], s=1)
    #ax.plot(rrr,bbb, c=[0.5]*4,linewidth=0.3)
    ax.set(xscale='log',yscale='log',xlabel='rho',ylabel='B')
    #ax.set(xlim=ext_d.minmax,ylim=ext_b.minmax)
    fig.savefig('%s/%s'%(plot_dir,fname))


if 'brho_rinf' not in dir():
    brho_rinf={}
if 'brho_r1' not in dir():
    brho_r1={}
    brho_rmax={}
if 'brho_rsmart_1' not in dir():
    brho_rsmart_1={}
if 'brho_rsmart_2' not in dir():
    brho_rsmart_2={}
if 1:
    for sim in sim_list:
        if sim in brho_rsmart_2:
            continue
        mon = monster.closet[sim]
        core_list =  mon.this_looper.core_by_mode['A']
        #this = meanie(mon)
        #this.run(core_list,frames='short',sphere_type='rsmart_1')
        #brho_rsmart_1[sim]=this
        this = meanie(mon)
        this.run(core_list,frames='movie',sphere_type='rsmart_2')
        brho_rsmart_2[sim]=this
        #this = meanie(mon)
        #this.run(core_list,frames='short',sphere_type='rinf')
        #brho_rinf[sim]=this
        #this = meanie(mon)
        #this.run(core_list, frames='short',sphere_type='r1')
        #brho_r1[sim]=this
        #that = meanie(mon)
        #that.run(core_list, frames='short',sphere_type='rmax')
        #brho_rmax[sim]=that

if 1:
    for sim in sim_list:
        ploot_bckgrnd(brho_rsmart_2[sim],fname='brho_rs2_tracks_%s.png'%sim)

if 0:
    #ploot(brho_rsmart_1[sim],fname='b_rho_%s_rsmart_1'%sim)
    #ploot(brho_rsmart_2[sim],fname='b_rho_%s_rsmart_2'%sim)
    ploot(brho_r1[sim],fname='b_rho_%s_r1'%sim)
    ploot(brho_rmax[sim],fname='b_rho_%s_rmax'%sim)
    ploot(brho_rinf[sim],fname='b_rho_%s_rinf'%sim)
    for sim in sim_list:
        ploot2(brho_r1[sim], brho_rmax[sim],vert='max',horz='1')
        ploot2(brho_r1[sim], brho_rinf[sim],vert='inf',horz='1')
        ploot2(brho_rmax[sim], brho_rinf[sim],vert='inf',horz='max')
        #ploot2(brho_rinf[sim], brho_rsmart[sim],vert='smart',horz='inf')
        #ploot2(brho_r1[sim], brho_rsmart[sim],vert='smart',horz='1')
        #ploot2(brho_rsmart_1[sim], brho_rsmart_2[sim],vert='smart_2',horz='smart_1')

