
from starter2 import *
import xtra_energy

import track_loader as TL

import pcolormesh_helper as pch
reload(pch)
def just_r_pseudo(mon,core_list,fracs):

    cores=[]
    rpseudo=[]
    for core_id in core_list:
        print('r pseudo',core_id)
        if core_id in [37,180]:
            print('core %d sucks'%core_id)
            continue
        frames = mon.frames_from_tsung(core_id,fracs)
        for nframe,frame in enumerate(frames[-1:]):
            proj_ax=0
            sph_small = mon.get_sphere(core_id,frame,'rinf')
            rho = np.abs(sph_small[YT_density])
            rfield = sph_small['radius']
            ok = rho > 1e4
            if (ok).any():
                rmax = rfield[ok].v.max()*colors.length_units_au
                cores.append(core_id)
                rpseudo.append(rmax)
                print(rmax)
    fig,ax=plt.subplots(1,1)
    print(rpseudo)
    import equal_probability_binner as epb
    #ax.hist(rpseudo)
    epb.equal_prob(nar(rpseudo),16,ax)
    fig.savefig('plots_to_sort/rpseudo')
    return nar(cores), nar(rpseudo)

def small_sphere(mon,core_list,fracs):

    for core_id in core_list:
        if core_id in [37,180]:
            print('core %d sucks'%core_id)
            continue
        frames = mon.frames_from_tsung(core_id,fracs)
        print('Sphere',core_id,frames)
        nrow = 3
        ncol = len(fracs)+2
        fig,ax=plt.subplots( nrow,ncol, figsize=(8,6))
        for nframe,frame in enumerate(frames):
            proj_ax=0
            sph8 = mon.get_sphere(core_id,frame,'r8')
            sph_small = mon.get_sphere(core_id,frame,'rinf')
            import r_inflection
            reload(r_inflection)
            rinf = mon.get_r_inflection(core_id,frame)
            #rinf,ue,rcen, due=r_inflection.RKEEP(sph8, return_ue=True)
            #rinf,ue,rcen, due=r_inflection.RKEEP(sph8, return_ue=True)
            ds   = mon.get_ds(frame)

            if 1:
                proj = ds.proj('density',proj_ax,center=sph8.center,data_source=sph8)
                dx = 1/2048
                Radius = 8./128
                Radius = ds.arr(Radius,'code_length')
                C = list(sph8.center)
                L = list(sph8.center-Radius)
                R = list(sph8.center+Radius)
                L.pop(proj_ax); R.pop(proj_ax);C.pop(proj_ax)
                ext=[L[0],R[0],L[1],R[1]]
                nz = int(2*Radius/dx)
                arr = proj.to_frb(2*Radius, nz)['density']
                norm = mpl.colors.LogNorm(vmin=arr[arr>0].min(),vmax=arr.max())
                circle = mpl.patches.Circle(C,radius=rinf, ec='red',fc='None')
                ax[0][nframe].imshow(arr,norm=norm,cmap='gray', extent=ext)
                ax[0][nframe].add_artist(circle)

            if 0:
                rho = np.abs(sph8[YT_grav_energy_2])
                rfield = sph8['radius']
                pch.simple_phase(rfield,rho,log=True,ax=ax[1][nframe])
                ax[1][nframe].set(xscale='log',yscale='log')
                ax[1][nframe].axvline(rinf)
                #axx = ax[1][nframe].twinx()
                #ax[1][nframe].plot( rcen, ue,c='r')
                #axx.plot( rcen[1:],due)
            if 1:
                rho = np.abs(sph_small[YT_density])
                rfield = sph_small['radius']
                pch.simple_phase(rfield,rho,log=True,ax=ax[2][nframe])
                ax[2][nframe].set(xscale='log',yscale='log', xlim=[1e-4,1e-1])
                ax[2][nframe].axvline(rinf, c='r')
                ok = rho > 1e4
                if (ok).any():
                    print('word')
                    rmax = rfield[ok].max()
                    ok2 = rfield == rmax
                    rho_ok = rho[ok2]
                    #print(rho_ok)
                    #print(rfield[ok2])
                    if rmax > 1e4/colors.length_units_au and frame==125:
                        pdb.set_trace()
                    #print(rmax*colors.length_units_au)
                    ax[2][nframe].axvline(rmax, c='g')


        fig.tight_layout()
        fig.savefig('plots_to_sort/rsmall_%s_c%04d.png'%(mon.name,core_id))


sims = ['u501']
TL.load_tracks(sims)
import monster
reload(monster)
monster.load(sims)

for sim in sims:
    L = TL.tracks[sim]
    core_list=L.core_by_mode['A']
    #core_list=[67]
    fracs_of_tsung = [0.0, 0.5]
    #small_sphere(monster.closet[sim],core_list,fracs_of_tsung)
    cores,rpseudo=just_r_pseudo(monster.closet[sim],core_list,fracs_of_tsung)




