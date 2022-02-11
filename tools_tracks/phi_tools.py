



from starter2 import *
import matplotlib.image as mpimg

from scipy.optimize import curve_fit
import data_locations as dl
from collections import defaultdict
import looper2
import heat_map
reload(heat_map)

G = 1620/(np.pi*4)
tff = np.sqrt(3*np.pi/32/G)
class phi_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.name = self.this_looper.sim_name
        self.cores_used=[]
        self.phi = defaultdict(list)
        self.Lsubsonic=defaultdict(list)
        self.Lsubvirial=defaultdict(list)
        self.Lmax = defaultdict(list)

        self.Msubsonic=defaultdict(list)
        self.Msubvirial=defaultdict(list)

        #ratios
        self.Lssx = defaultdict(list)
        self.Lsvx = defaultdict(list)
        self.Lsssv = defaultdict(list)

    def run(self,core_list=None, do_plots=None):

        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        self.times=tsorted

        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id)
            tmap=rainbow_map(ms.ntimes)
            if ms.nparticles < 10:
                continue
            print('SubLength ', core_id)
            self.cores_used.append(core_id)
            n0=0
            rmin = 1./2048
            fig,ax=plt.subplots(1,2, figsize=(8,4))
            rmap = rainbow_map(ms.ntimes)
            phi_all = thtr.c([core_id],'PotentialField')


            if 0:
                fig,ax=plt.subplots(1,1)
                for nt, time in enumerate(thtr.times):
                    density = thtr.c([core_id],'density')[:,nt]
                    phi = thtr.c([core_id],'PotentialField')[:,nt]
                    cell_volume = thtr.c([core_id],'cell_volume')[:,nt]

                    vals , bins = np.histogram(phi,weights=cell_volume,density=True)
                    bc = 0.5*(bins[1:]+bins[:-1])

                    ax.plot( bc,vals, c = tmap(nt))
                fig.savefig('plots_to_sort/%s_phi_c%04d'%(self.this_looper.sim_name, core_id))

            if 1:
                velocity_method = 'vtc'
                if do_plots:
                    fig,axbox=plt.subplots(2,2, figsize=(12,12))
                    ax=axbox.flatten()
                ek_ext = extents()
                eg_ext = extents()
                abs_max=extents()
                vir_ext=extents()
                for nt, time in enumerate(thtr.times):
                    ms.get_central_velocity(core_id,nt)

                    density = thtr.c([core_id],'density')[:,nt]
                    phi = thtr.c([core_id],'PotentialField')[:,nt]
                    cell_volume = thtr.c([core_id],'cell_volume')[:,nt]
                    if velocity_method == 'vrm':
                        vr = ms.vr_rel[:,nt]  #the radial one, works ok
                    if velocity_method == 'vtm':
                        vr = ms.rel_vmag[:,nt]  #testing
                    if velocity_method == 'vtc':
                        vr = ms.cen_vmag[:,nt]  #testing
                    if velocity_method == 'vrc':
                        vr = ms.rc_vmag[:,nt]  #testing
                    if velocity_method == 'vtcd':
                        vr = ms.rcd_vmag[:,nt]  #testing
                    this_r = ms.r[:,nt]
                    this_r[ this_r < rmin] = rmin
                    asort = np.argsort(this_r)

                    EG_interior = G*density*phi*cell_volume
                    EK_interior = G*density*vr**2*cell_volume
                    vir=EK_interior/np.abs(EG_interior)
                    if do_plots:
                        ax[0].scatter( this_r, -EG_interior, c=[tmap(nt)]*len(this_r))
                        ax[1].scatter( this_r, EK_interior, c=[tmap(nt)]*len(this_r))
                        ax[2].scatter( this_r, vir, c=[tmap(nt)]*len(this_r))
                    ek_ext(EK_interior)
                    eg_ext(EG_interior)
                    abs_max(np.abs(EK_interior))
                    abs_max(np.abs(EG_interior))
                    vir_ext(vir)

                    vrs = vr[asort]
                    dv = cell_volume[asort]
                    sigma_vr2 = np.cumsum(vrs**2*dv)/np.cumsum(dv)
                    vir_sorted = vir[asort]
                    r_sorted = this_r[asort]
                    Lsubsonic = 0
                    moose = (density*cell_volume)[asort] #its not the mass.
                    mooses=moose.sum()
                    if (sigma_vr2 < 1).any():
                        subsonic_index = max(np.where( sigma_vr2[asort] <= 1)[0])
                        Lsubsonic = r_sorted[subsonic_index]
                        Msubsonic = moose[subsonic_index]/mooses
                    if (vir < 1).any():
                        subvir_index = max(np.where( vir[asort] <= 1)[0])
                        Lsubvirial = r_sorted[subvir_index]
                        Msubvirial = moose[subsonic_index]/mooses
                    

                    self.Msubsonic[core_id].append(Msubsonic)
                    self.Msubvirial[core_id].append(Msubvirial)
                    Lmax= r_sorted.max()
                    self.Lsubsonic[core_id].append(Lsubsonic)
                    self.Lsubvirial[core_id].append(Lsubvirial)
                    self.Lmax[core_id].append(Lmax)

                    self.Lssx[core_id].append( Lsubsonic/Lmax)
                    self.Lsvx[core_id].append( Lsubvirial/Lmax)
                    r3= Lsubsonic/Lsubvirial
                    self.Lsssv[core_id].append(r3)

                    #vals,bins=np.histogram(vir)
                    #bincen=0.5*(bins[1:]+bins[:-1])
                    #ax[3].plot( bincen, vals, c=tmap(nt))



                    #vals , bins = np.histogram(phi,weights=cell_volume,density=True)
                    #bc = 0.5*(bins[1:]+bins[:-1])

                    #ax.plot( bc,vals, c = tmap(nt))
                #fig.savefig('plots_to_sort/%s_phi_c%04d'%(self.this_looper.sim_name, core_id))
                if do_plots:
                    axbonk(ax[0],xscale='log',xlabel='r',ylabel='EG')
                    axbonk(ax[1],xscale='log',xlabel='r',ylabel='EK')
                    axbonk(ax[2],xscale='log',xlabel='r',ylabel='VIR')
                    ax[0].set_yscale('symlog',linthresh=0.002)
                    ax[0].set_ylim(-1e-3,abs_max.minmax[1])
                    ax[1].set_yscale('symlog',linthresh=0.002)
                    ax[1].set_ylim(-1e-3,abs_max.minmax[1])
                    max_vir=max(vir_ext.minmax)
                    ax[2].set_yscale('symlog', linthresh=max_vir/100)
                    ax[2].set_ylim(-max_vir,max_vir)

                    #axbonk(ax[3],xlabel='Virial',ylabel='N',xscale='log',yscale='log')
                    ax[3].plot( thtr.times/tff, self.Lsubsonic[core_id], c='r')
                    ax[3].plot( thtr.times/tff, self.Lsubvirial[core_id], c='g')
                    ax[3].plot( thtr.times/tff, self.Lmax[core_id], c='k')
                    axbonk(ax[3], ylabel='L',xlabel=r'$t/t_{ff}$')
                    #axbonk(ax[3],xlabel='EK',ylabel='EG')
                    #ax[3].set_xscale('symlog',linthresh=0.0002)
                    #ax[3].set_yscale('symlog',linthresh=0.0002)
                    fig.savefig('plots_to_sort/%s_EG_c%04d'%(self.this_looper.sim_name, core_id))



import three_loopers_tenfour as TL4

if 'Ptools' not in dir():
    Ptools={}
if 0:
    sim_list = ['u401']
    for this_simname in sim_list:
        if this_simname not in Ptools:
            Ptools[this_simname] = phi_tool(TL4.loops[this_simname])
            Ptools[this_simname].run(do_plots=False)

sim_list=['u502']
if 'these_loops' not in dir():
    these_loops={}

if 'u502' not in Ptools:
    u500={'u501':"/data/cb1/Projects/P19_CoreSimulations/CoreSets/u500/u501_all_frame_all_prim.h5",
          'u502':"/data/cb1/Projects/P19_CoreSimulations/CoreSets/u500/u502_all_frame_all_prim.h5"
    for sim in sim_list:
        print('LOAD LOOP u501')
        these_loops['u501']=looper2.load_looper(u501)
        Ptools['u501'] = phi_tool( these_loops['u501'])
        Ptools['u501'].run(do_plots=False)

for this_simname in sim_list:
    fig,ax=plt.subplots(2,3,figsize=(12,8))
    axlist=ax.flatten()
    qmat = heat_map.plot_heat(tool=Ptools[this_simname], quan_dict=Ptools[this_simname].Lmax, ax=axlist[0])
    axbonk(axlist[0],xlabel='t/tff',ylabel='Lmax', ylim=[0,0.65])
    qmat = heat_map.plot_heat(tool=Ptools[this_simname], quan_dict=Ptools[this_simname].Lsubvirial, ax=axlist[1])
    axbonk(axlist[1],xlabel='t/tff',ylabel='Lsubvirial', ylim=[0,0.65])
    qmat = heat_map.plot_heat(tool=Ptools[this_simname], quan_dict=Ptools[this_simname].Lsubsonic, ax=axlist[2])
    axbonk(axlist[2],xlabel='t/tff',ylabel='Lsubsonic', ylim=[0,0.65])

    LsLvExt = extents()
    for core_id in Ptools[this_simname].Lsssv:
        arr=nar(Ptools[this_simname].Lsssv[core_id])
        LsLvExt( arr[arr>0])
    LsLxExt = extents()
    for core_id in Ptools[this_simname].Lssx:
        arr=nar(Ptools[this_simname].Lssx[core_id])
        LsLxExt( arr[arr>0])
    LvLxExt = extents()
    for core_id in Ptools[this_simname].Lsvx:
        arr=nar(Ptools[this_simname].Lsvx[core_id])
        LvLxExt( arr[arr>0])


    bins = np.geomspace( 1e-3,40,64)
    stuff_lsssv = heat_map.plot_heat(tool=Ptools[this_simname], quan_dict=Ptools[this_simname].Lsssv, ax=axlist[3],bins=bins)
    axbonk(axlist[3],xlabel='t/tff',ylabel='Lss/Lsv')
    axlist[3].set_yscale('symlog',linthresh=0.1)
    axlist[3].set_ylim([0,40])
    times=Ptools[this_simname].this_looper.tr.times/tff
    axlist[3].plot(times, times*0+1,c='k')
    bins = np.geomspace( 1e-3,1.01,64)
    stuff_lsvx = heat_map.plot_heat(tool=Ptools[this_simname], quan_dict=Ptools[this_simname].Lsvx, ax=axlist[4],bins=bins)
    axbonk(axlist[4],xlabel='t/tff',ylabel='Lsubvirial/Lmax')
    axlist[4].set_yscale('symlog',linthresh=0.1)
    axlist[4].plot(times, times*0+1/2,c='k')

    stuff_lssx = heat_map.plot_heat(tool=Ptools[this_simname], quan_dict=Ptools[this_simname].Lssx, ax=axlist[5],bins=bins)
    axlist[5].plot(times, times*0+1/2,c='k')
    axbonk(axlist[5],xlabel='t/tff',ylabel='Lsubsonic/Lmax')
    axlist[5].set_yscale('symlog',linthresh=0.1)

    fig.savefig('plots_to_sort/%s_LxLsLv.png'%this_simname)
    plt.close('all')

    fig,ax=plt.subplots(1,2)
    fracs = np.linspace(0,1,11)
    frac_map = rainbow_map(fracs.size)
    sets_lsvx=[]
    TheY_lsvx = stuff_lsvx['TheY']
    TheH_lsvx = stuff_lsvx['hist']
    dv_lsvx = stuff_lsvx['dv']
    dv_lsvx.shape = 1,dv_lsvx.size

    sets_lssx=[]
    TheY_lssx = stuff_lssx['TheY']
    TheH_lssx = stuff_lssx['hist']
    dv_lssx = stuff_lssx['dv']
    dv_lssx.shape = 1,dv_lssx.size
    for nf,fr in enumerate(fracs[:-1]):
        mask = ( (TheY_lsvx > fracs[nf])*(TheY_lsvx <= fracs[nf+1]) ).astype('int')
        sets_lsvx.append( (TheH_lsvx*dv_lsvx*mask).sum(axis=1))

        mask = ( (TheY_lssx > fracs[nf])*(TheY_lssx <= fracs[nf+1]) ).astype('int')
        sets_lssx.append( (TheH_lssx*dv_lssx*mask).sum(axis=1))
    ax[0].stackplot(times, *sets_lsvx)#, c=frac_map(nf))
    ax[1].stackplot(times, *sets_lssx)#, c=frac_map(nf))
    fig.savefig('plots_to_sort/Lsub_stack.png')








