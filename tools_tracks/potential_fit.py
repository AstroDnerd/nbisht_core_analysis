from starter2 import *
import xtra_energy
import three_loopers_u500 as TL
import pcolormesh_helper as pch 
sim_list=['u501']

def plot_phi(this_looper,core_list=None):
    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    frame = this_looper.target_frame
    for core_id in core_list:
        print('Potential %s %d'%(this_looper.sim_name,core_id))
        save_prefix='plots_to_sort/%s_n%04d_c%04d'%(this_looper.sim_name,frame,core_id)
        ds = this_looper.load(frame)
        xtra_energy.add_energies(ds)
        ms = trackage.mini_scrubber(this_looper.tr,core_id)
        c = nar([ms.mean_x[-1], ms.mean_y[-1],ms.mean_z[-1]])
        
        rsph = ds.arr(8.0/128,'code_length')
        #rsph = ds.arr(2/128,'code_length')
        sp = ds.sphere(c,rsph)

        if 0:
            proj = ds.proj('grav_energy',0, center=c,data_source=sp)
            pw = proj.to_pw(center=c)
            pw.zoom(32)
            pw.save(save_prefix)

        GE = np.abs(sp['grav_energy'])
        dv = np.abs(sp['cell_volume'])
        RR = sp['radius']
        r_sphere=sp['radius']
        DD = sp['density']

        ok = RR < 0.01

        fit = np.polyfit( np.log10(RR[ok]), np.log10(DD[ok]), 1)
        M = sp['cell_mass'].sum()
        G = 1620/(4*np.pi)
        coe = -1/(8*np.pi*G)
        power=2*fit[0]+2
        phi_del_squ_analy = (coe*G**2*M**2*r_sphere**(power))/max(r_sphere)**(2*fit[0]+6)


        fig,ax=plt.subplots(1,2)
        ax0=ax[0];ax1=ax[1]
        if 1:
            rbins = np.geomspace( r_sphere[r_sphere>0].min(), r_sphere.max(),64)
            gbins = np.geomspace( GE[GE>0].min(), GE.max(),64)
            hist, xb, yb = np.histogram2d( r_sphere, GE, bins=[rbins,gbins],weights=dv)
            pch.helper(hist,xb,yb,ax=ax0,transpose=False)
        if 1:
            rbins = np.geomspace( r_sphere[r_sphere>0].min(), r_sphere.max(),64)
            dbins = np.geomspace( DD[DD>0].min(), DD.max(),64)
            hist, xb, yb = np.histogram2d( r_sphere, DD, bins=[rbins,dbins],weights=dv)
            pch.helper(hist,xb,yb,ax=ax1,transpose=False)
        if 0:
            ax.scatter( sp['radius'],-sp['grav_energy'],c='r',s=0.1)
            #ax.set_xlim(sp['radius'].min(),sp['radius'].max())

        if 1:
            axt = ax0.twinx()
            sort = np.argsort(sp['radius'])
            rsort = sp['radius'][sort]
            GEdv = (GE*RR)[sort]
            tots = np.cumsum(GEdv)
            axt.plot( rsort, tots)


        ax0.plot(  RR, np.abs(phi_del_squ_analy),c='k')
        ax1.plot( RR[ok], 10**( fit[0]*np.log10(RR[ok])+fit[1]),c='r')
        axbonk(ax0,xscale='log',yscale='log',xlabel='r',ylabel='abs(grav_eng)')
        axbonk(ax1,xscale='log',yscale='log',xlabel='r',ylabel='density')
        fig.savefig(save_prefix+"grav_eng")
        continue
        ax.clear()
        bins = np.geomspace( GE[GE>0].min(), GE.max(),64)
        ax.hist( np.abs(sp['grav_energy']), bins=bins, histtype='step')
        axbonk(ax,xlabel='Grad Phi',ylabel='N',xscale='log',yscale='log')
        fig.savefig(save_prefix+"grad_phi_hist")


        min_ge=np.argmin(GE)
        x,y,z = sp['x'][min_ge], sp['y'][min_ge], sp['z'][min_ge]
        c = ds.arr([x,y,z],'code_length')
        SSS = yt.SlicePlot(ds,'x','grav_energy',center=c)
        SSS.annotate_sphere(center=c,radius=rsph)
        SSS.save(save_prefix)
        SSS = yt.SlicePlot(ds,'x','kinetic_energy',center=c)
        #SSS.annotate_sphere(sp)
        SSS.save(save_prefix)

all_cores=np.unique( TL.loops['u501'].tr.core_ids)
plot_phi( TL.loops['u501'])#,core_list=[323]+list(all_cores)[0:3])


