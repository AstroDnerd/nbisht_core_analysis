from starter2 import *
import xtra_energy
reload(xtra_energy)
import three_loopers_u500 as TL
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import pcolormesh_helper as pch
import colors
import movie_frames 
from collections import defaultdict
from scipy.ndimage import gaussian_filter
import other_scrubber
reload(other_scrubber)
G = colors.G
class vel_color():
    def __init__(self):
        self.alpha_time=defaultdict(list)
        self.cores_used=[]
    def run(self,this_looper,core_list=None, do_plots=True, r_inflection=None, frame_list=None, tsing=None):

        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        thtr=this_looper.tr
        if frame_list is None and False:
            mask = movie_frames.quantized_mask(this_looper).flatten()
            ok = np.zeros_like(mask)
            ok[::10] = mask[::10]
            mask=ok ;print('kludge mask')
            times=thtr.times[mask]+0 #the zero makes a copy
            self.times=times+0
            times.shape=times.size,1
            times=times/colors.tff
            frame_list=thtr.frames[mask]
            rm = rainbow_map(len(frame_list))




        fig5,ax5=plt.subplots(3,2)
        y_ext=extents()
        for core_id in core_list:
            if 1:
                just_tsing = np.zeros_like(thtr.times, dtype='bool')
                index_sing=np.argmin( np.abs( thtr.times/colors.tff-tsing.tsing_core[core_id]))
                index_end=np.argmin( np.abs( thtr.times/colors.tff-tsing.tend_core[core_id]))
                #print(tsing.tsing_core[core_id])
                #print(tsing.tend_core[core_id])
                just_tsing[index_sing]=True
                #just_tsing[index_end]=True
                times = thtr.times[just_tsing]+0
                times.shape=times.size,1
                frame_list=thtr.frames[just_tsing]
                rm = rainbow_map(len(frame_list))
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            self.cores_used.append(core_id)
            
            r_bins = np.geomspace( 2e-4, 32/128, 32)
            r_cen = 0.5*(r_bins[1:]+r_bins[:-1])
            RV = np.zeros( [len(r_bins)-1, len(frame_list)])
            frame_rmap=rainbow_map(len(frame_list))
            fig6,ax6=plt.subplots(1,3,figsize=(12,8))
            ax1=ax6[0]; ax2=ax6[1]
            got_tsing=False
            got_tend=False
            for nframe,frame in enumerate(frame_list):
                print("profile on %s c%04d n%04d"%(this_looper.sim_name, core_id, frame))
                ds = this_looper.load(frame)
                xtra_energy.add_energies(ds)
                xtra_energy.add_gdotgradrho(ds)
                nf = np.where( this_looper.tr.frames == frame)[0][0]
                time = times[nframe,0]

                c = nar([ms.mean_xc[nf], ms.mean_yc[nf],ms.mean_zc[nf]])
                msR = ms.rc
                msR[ msR<1/2048]=1/2048
                
                MaxRadius=msR[:,nf].max()
                Radius = max([8.0/128, MaxRadius])
                rsph = ds.arr(Radius,'code_length')
                sp = ds.sphere(c,rsph)
                dv = sp[YT_cell_volume]
                RR = sp['radius']
                DD = sp[YT_density]
                ORDER = np.argsort( RR)
                RR_cuml = RR[ORDER]
                rho_srt = DD[ORDER]
                V_cuml = np.cumsum( dv[ORDER])
                M_cuml = np.cumsum( DD[ORDER]*dv[ORDER])
                dv_sort = dv[ORDER]

                r_bins=np.geomspace(RR_cuml[10],RR.max(),64)
                r_cen=0.5*(r_bins[1:]-r_bins[:-1])
                digitized = np.digitize( RR_cuml, r_bins)
                def digit(arr):
                    output  =nar([ arr[ digitized == i].mean() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])
                    return output


                rho_bin = digit(rho_srt)
                ok = ~np.isnan(rho_bin)*(rho_bin>0)
                rrr = digit(RR_cuml)


                pfit_mass = np.polyfit( np.log10(RR_cuml[1:]), np.log10(M_cuml[1:]),1)
                pfit_rho  = np.polyfit( np.log10(RR_cuml[1:]), np.log10(rho_srt[1:]),1)
                pfit_dig = np.polyfit( np.log10(rrr[ok]), np.log10(rho_bin[ok]),1)
                pfit_m2  = np.polyfit( np.log10(RR_cuml[1:]), np.log10((M_cuml/V_cuml))[1:],1)


                ax1.plot( RR_cuml, rho_srt, c=[0.5]*4)
                ax1.plot( rrr[ok], rho_bin[ok],c='k')
                ax1.plot( RR_cuml, M_cuml/V_cuml,'b')

                ax1.plot( rrr[ok], 10**(np.log10(rrr[ok])*pfit_dig[0]+pfit_dig[1]),'k:',label="dig %0.2f"%pfit_dig[0])
                ax1.plot( rrr[ok], 10**(np.log10(rrr[ok])*pfit_rho[0]+pfit_rho[1]),'r', label="fit %0.2f"%pfit_rho[0])
                ax1.plot( rrr[ok], 10**(np.log10(rrr[ok])*pfit_m2[0]+pfit_m2[1]),'b:', label="m2 %0.2f"%pfit_m2[0])
                ax1.legend(loc=0)

                ax2.plot( RR_cuml, M_cuml,c=[0.5]*3)
                ax2.plot( RR_cuml, 10**(pfit_mass[0]*np.log10(RR_cuml)+pfit_mass[1]),"k--",label='fit %0.2f'%( pfit_mass[0]))
                MMM = pfit_rho[0]+3
                BBB = pfit_mass[1]
                ax2.plot( RR_cuml, 10**(MMM*np.log10(RR_cuml)+BBB),label='a+3: %0.2f'%(MMM), c='r')
                MMM = pfit_dig[0]+3
                BBB = pfit_mass[1]
                ax2.plot( RR_cuml, 10**(MMM*np.log10(RR_cuml)+BBB),label='dig+3: %0.2f'%(MMM), c='k')
                ax2.set(title='M')
                ax2.legend(loc=0)

                d2_sort = DD[ORDER]**2
                d2_cuml = np.cumsum( DD[ORDER]**2*dv[ORDER])#/V_cuml
                pfit_d2 = np.polyfit( np.log10(RR_cuml[1:]), np.log10( d2_cuml[1:]),1)
                ax6[2].plot(RR_cuml, d2_sort, c=[0.5]*4)
                ax6[2].plot(RR_cuml[5:], d2_cuml[5:], c='k')
                ax6[2].plot( RR_cuml, 10**(np.log10(RR_cuml)*pfit_d2[0]+pfit_d2[1]),'k:', label="d2 %0.2f"%pfit_d2[0])

                ax6[2].plot( RR_cuml, (M_cuml/V_cuml)**2, c='g')
                ax6[2].legend(loc=0)
                ax6[2].set(yscale='log',xscale='log')


#ax6[1][0].plot(RR_cuml, rho_toy*4*np.pi*RR_cuml**3)
                for aaa in ax6.flatten():
                    aaa.set(yscale='log',xscale='log')

                fig6.savefig('plots_to_sort/rho_squared_%s_c%04d_n%04d'%(this_looper.sim_name, core_id, frame))
import three_loopers_u500 as TL
sim_list=['u501','u502','u503']
sim_list=['u502']

import tsing
reload(tsing)
if 'tsing_tool' not in dir() or True:
    tsing_tool={}
    for ns,sim in enumerate(sim_list):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()

if 'RV' not in dir() or True:
    for sim in sim_list:
        all_cores=np.unique( TL.loops[sim].tr.core_ids)
        core_list=[323]
        core_list=[25]
        core_list=[74]
        #core_list = TL.loops[sim].core_by_mode['Alone']
        #core_list=core_list[10:]
        #core_list=[114]
        TF=TL.loops[sim].target_frame
        all_frames = TL.loops[sim].tr.frames
        nframes = all_frames.size
        frame_list = all_frames #[0, TF, all_frames[int(nframes/2)]]
        #frame_list = [100]

        RV=vel_color()
        RV.run(TL.loops[sim],core_list=core_list, do_plots=True, frame_list=None, tsing=tsing_tool[sim])#, r_inflection=anne.inflection[sim])
