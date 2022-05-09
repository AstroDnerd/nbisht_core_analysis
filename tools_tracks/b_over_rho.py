from starter2 import *
import xtra_energy
reload(xtra_energy)

import three_loopers_u500 as TL
if 1:
    for sim in ['u501']:
        this_looper=TL.loops[sim]
        fig,ax=plt.subplots(1,1)
        all_frames=this_looper.tr.frames[::10]
        frames=[0,50,all_frames[-1]]
        core_list=[323]
        for frame in frames:
            ds = this_looper.load(frame)
            xtra_energy.add_b_over_rho(ds)
            for core_id in core_list:
                ms = trackage.mini_scrubber( this_looper.tr, core_id)
                all_c = np.stack([ms.mean_x,ms.mean_y,ms.mean_z])
                proj = ds.proj(YT_bx_over_rho,0,center=c)
                pw=proj.to_pw()
                #pw.zoom(128)
                print(pw.save('plots_to_sort/%s_n%04d'%(sim,frame)))
if 0:
    for sim in ['u501']:
        this_looper=TL.loops[sim]
        fig,ax=plt.subplots(1,1)
        rm = rainbow_map(this_looper.tr.frames.size)
        all_frames=this_looper.tr.frames[::10]
        frames=[0,50,all_frames[-1]]
        for frame in frames:
            ds = this_looper.load(frame)
            xtra_energy.add_b_over_rho(ds)
            ad = ds.all_data()
            dxmin=ds.index.get_smallest_dx()
            bx_over_rho = ad[YT_bx_over_rho]/np.abs(ad['velocity_magnitude']/dxmin)
            #by_over_rho = ad[YT_by_over_rho]/np.abs(ad['velocity_magnitude']/dxmin)
            #bz_over_rho = ad[YT_bz_over_rho]/np.abs(ad['velocity_magnitude']/dxmin)
            ax.hist( bx_over_rho.v, label='bx',histtype='step',color= rm(frame))
            #ax.hist( by_over_rho.v, label='by',histtype='step',color= rm(frame))
            #ax.hist( bz_over_rho.v, label='bz',histtype='step',color= rm(frame))
        ax.set_yscale('log')
        fig.savefig('plots_to_sort/bi_over_rho_%s'%sim)


    if 0:
        proj = ds.proj(YT_bx_over_rho,0)
        pw=proj.to_pw()
        pw.save('plots_to_sort/%s_n%04d_all'%(sim,frame))
