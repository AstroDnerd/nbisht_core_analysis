
from starter2 import *
import track_loader as TL
import xtra_energy
import pcolormesh_helper as pch

def extract(arr, sl, mask):
    return arr[sl].transpose()[mask,:]
class get_rb():
    def __init__(self,sphere=None,thtr=None,core_id=None, sl=None, mask=None):
        if thtr is not None:
            bx   =extract(thtr.c([core_id],'magnetic_field_x'),sl,mask)
            by   =extract(thtr.c([core_id],'magnetic_field_y'),sl,mask)
            bz   =extract(thtr.c([core_id],'magnetic_field_z'),sl,mask)
            vx   =extract(thtr.c([core_id],'velocity_x'),sl,mask)
            vy   =extract(thtr.c([core_id],'velocity_y'),sl,mask)
            vz   =extract(thtr.c([core_id],'velocity_z'),sl,mask)
            dxvx =extract(thtr.c([core_id],'dxvx'),sl,mask)
            dxvy =extract(thtr.c([core_id],'dxvy'),sl,mask)
            dxvz =extract(thtr.c([core_id],'dxvz'),sl,mask)
            dyvx =extract(thtr.c([core_id],'dyvx'),sl,mask)
            dyvy =extract(thtr.c([core_id],'dyvy'),sl,mask)
            dyvz =extract(thtr.c([core_id],'dyvz'),sl,mask)
            dzvx =extract(thtr.c([core_id],'dzvx'),sl,mask)
            dzvy =extract(thtr.c([core_id],'dzvy'),sl,mask)
            dzvz =extract(thtr.c([core_id],'dzvz'),sl,mask)
            B2   =extract(thtr.c([core_id],'magnetic_field_strength'),sl,mask)**2
            divv =extract(thtr.c([core_id],'velocity_divergence'),sl,mask)
        if sphere is not None:
            bx   =sphere['magnetic_field_x']
            by   =sphere['magnetic_field_y']
            bz   =sphere['magnetic_field_z']
            vx   =sphere['velocity_x']
            vy   =sphere['velocity_y']
            vz   =sphere['velocity_z']
            dxvx =sphere['dxvx']
            dxvy =sphere['dxvy']
            dxvz =sphere['dxvz']
            dyvx =sphere['dyvx']
            dyvy =sphere['dyvy']
            dyvz =sphere['dyvz']
            dzvx =sphere['dzvx']
            dzvy =sphere['dzvy']
            dzvz =sphere['dzvz']
            B2   =sphere['magnetic_field_strength']**2
            divv =sphere['velocity_divergence']
        Sx = bx*dxvx+by*dyvx+bz*dzvx
        Sy = bx*dxvy+by*dyvy+bz*dzvy
        Sz = bx*dxvz+by*dyvz+bz*dzvz
        Stretch= bx*Sx+by*Sy+bz*Sz
        RB=Stretch/(B2*divv)
        self.Sx=Sx
        self.Sy=Sy
        self.Sz=Sz
        self.B2=B2
        self.divv=divv
        self.RB=RB




class kappa_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.B = defaultdict(list)
        self.rho=defaultdict(list)

        self.meanRb = defaultdict(list)

        self.pB = defaultdict(list)
        self.prho=defaultdict(list)

    def run(self, core_list=None,frame_list=None, tsing=None):
        this_looper=self.this_looper
        thtr=this_looper.tr


        if core_list is None:
            core_list=np.unique(this_looper.tr.core_ids)
        if frame_list is None:
            frame_list=this_looper.tr.frames

        self.cores_used=copy.copy(core_list)
        self.times=[]

        first_light_mask={}
        scrubs={}
        for core_id in core_list:
            ms=trackage.mini_scrubber(thtr,core_id, do_magnetic=True)
            scrubs[core_id]=ms
            tsung = tsing.tend_core[core_id]
            nf = np.where(thtr.times/colors.tff==tsung)[0][0]
            ok = ms.density[:,nf]>1e4
            first_light_mask[core_id]=ok
        for iframe, frame in enumerate(frame_list):
            ds = this_looper.load(frame)
            xtra_energy.add_v_grad(ds)
            radius=1e-2
            radius = ds.arr(radius,'code_length')
            nframe=np.where(thtr.frames==frame)[0][0]
            self.times.append(thtr.times[nframe])

            los = 0
            XYZ = 'xyz'[los]
            BI = 'magnetic_field_'+XYZ
            BI = 'magnetic_field_strength'

            for core_id in core_list:
                print('yay',core_id, frame)
                ms = scrubs[core_id]
                p = np.array([ms.mean_x[nframe],ms.mean_y[nframe],ms.mean_z[nframe]])
                c = ds.arr(p,'code_length')
                r_ok= ms.r[first_light_mask[core_id],nframe]
                r = max([ 1./128, r_ok.max()])


                sph = ds.sphere(c,r)

                rrr = sph['radius']
                dv = sph['cell_volume']
                B  = sph[BI]
                rho = sph['density']

                self.B[core_id].append( (B*dv*rho).sum()/(dv*rho).sum())
                self.rho[core_id].append( (rho*dv).sum()/(dv).sum())

                mask = ms.compute_unique_mask(core_id, 1/2048,nframe)
                pB = np.sqrt(ms.b2[mask,nframe])
                prho= ms.density[mask,nframe]
                dv  = ms.cell_volume[mask,nframe]
                pBmean = (pB*prho*dv).sum()/(prho*dv).sum()
                prhomean=(prho*dv).sum()/dv.sum()


                self.pB[core_id].append(pBmean)
                self.prho[core_id].append(prhomean)

                rbstuff = get_rb(sph)
                Rb = rbstuff.RB
                meanRb= Rb[np.abs(Rb)<10].mean() 
                self.meanRb[core_id].append(meanRb)

                if 1:
                    fig,axes=plt.subplots(2,2)
                    ax0=axes[0][0];ax1=axes[0][1];ax2=axes[1][0];ax3=axes[1][1]

                    proj=ds.proj(YT_density,1,data_source=sph, center=c)
                    frb = proj.to_frb(2*r,512)
                    density=frb[YT_density]
                    ax0.imshow(np.log10(density))
                    field=frb[YT_magnetic_field_strength]
                    ax1.imshow(field)

                    density = sph[YT_density]
                    field = sph[YT_magnetic_field_strength]
                    rho_bins = np.geomspace(density[density>0].min(),density.max(),16)
                    field_bins = np.geomspace(field[field>0].min(),field.max(),16)
                    pch.simple_phase(density.flatten(),field.flatten(),bins=[rho_bins,field_bins],ax=ax2)
                    ok = (density>0)*(field>0)
                    pfit = np.polyfit(np.log10(density[ok].flatten()), np.log10(field[ok].flatten()),1)
                    ax2.plot(rho_bins, 10**(pfit[0]*np.log10(rho_bins)+pfit[1]))
                    kappa_rb=1-meanRb
                    print("kappa_rb",kappa_rb)
                    ax2.plot(rho_bins, 10**(kappa_rb*np.log10(rho_bins)))
                    ax2.set(xscale='log',yscale='log',title="%0.3f"%pfit[0], ylim=[field[field>0].min(),field.max()])


                    density = frb[YT_density]
                    field = frb[YT_magnetic_field_strength]
                    rho_bins = np.geomspace(density[density>0].min(),density.max(),16)
                    field_bins = np.geomspace(field[field>0].min(),field.max(),16)
                    pch.simple_phase(density.flatten(),field.flatten(),bins=[rho_bins,field_bins],ax=ax3)
                    ok = (density>0)*(field>0)
                    pfit = np.polyfit(np.log10(density[ok].flatten()), np.log10(field[ok].flatten()),1)
                    ax3.plot(rho_bins, 10**(pfit[0]*np.log10(rho_bins)+pfit[1]))
                    ax3.set(xscale='log',yscale='log',title="%0.3f"%pfit[0])
                    fig.savefig('%s/%s_spheres_c%04d_n%04d'%(plot_dir,this_loop.sim_name,core_id,frame))



sim_list = ['u902']
TL.load_tracks(sim_list)

import tsing
tsing_tool=tsing.get_tsing(TL.tracks)

if 'rrr' not in dir() or True:
    rrr = {}

if 1:
    for sim in sim_list:
        if sim not in rrr:
            this_loop=TL.loops[sim]
            rr = kappa_tool(this_loop)
            #frame_list=this_loop.tr.frames[::10]
            frame_list=this_loop.tr.frames[-1:]
            frame_list=this_loop.tr.frames[-8:-7]
            core_list = TL.loops[sim].core_by_mode['Alone'][:2]
            rr.run(core_list=core_list,frame_list=frame_list, tsing=tsing_tool[sim])
            rrr[sim]=rr

if 1:
    for sim in sim_list:
        rr = rrr[sim]
        fig,axes=plt.subplots(1,2)
        ax0=axes[0];ax1=axes[1]
        for core_id in rr.cores_used:
            tsung = tsing_tool[sim].tend_core[core_id]
            print(tsung,core_id)
            ax0.plot( rr.times/tsung/colors.tff, np.log(rr.B[core_id]/rr.B[core_id][0]),c='r')
            ax0.plot( rr.times/tsung/colors.tff, np.log(rr.rho[core_id]/rr.rho[core_id][0]),c='k')
            kappa = np.log(rr.B[core_id]/rr.B[core_id][0])/np.log(rr.rho[core_id]/rr.rho[core_id][0])
            ax1.plot( rr.times, kappa, c='k')
            Rb = nar(rr.meanRb[core_id])
            ax1.plot(rr.times, 1-Rb, c='r')
        ax0.set(yscale='linear')
        ax1.set(ylim=[0,1])
        outname='%s/%s_kappa_rb'%(plot_dir,sim)
        fig.savefig(outname)
        print('save',outname)

if 0:

    fig,ax=plt.subplots(2,2)
    ax0=ax[0][0];ax1=ax[0][1]; ax2=ax[1][0];ax3=ax[1][1]
    
    the_x,the_y=np.log(rrr.rho), np.log(rrr.B)
    fit = np.polyfit( the_x,the_y,1)
    ax0.plot(the_x,fit[0]*the_x+fit[1])
    ax0.scatter(the_x,the_y)
    ax0.set_title("2d %0.2f"%fit[0])

    the_x,the_y=np.log(rrr.prho), np.log(rrr.pB)
    fit = np.polyfit( the_x,the_y,1)
    ax1.plot(the_x,fit[0]*the_x+fit[1])
    ax1.scatter(the_x,the_y)
    ax1.set_title("3d %0.2f"%fit[0])

    ax2.scatter( rrr.rho, rrr.prho)
    ax3.scatter( rrr.B, rrr.pB)
    ax2.set(xscale='log',yscale='log')
    ax3.set(xscale='log',yscale='log')



    fig.savefig('plots_to_sort/derp.png')
    plt.close('all')
