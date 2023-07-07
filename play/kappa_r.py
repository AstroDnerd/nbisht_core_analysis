from starter2 import *
import xtra_energy
import track_loader as TL
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import pcolormesh_helper as pch
import colors
import movie_frames 
import other_scrubber
reload(other_scrubber)
G = colors.G
plt.close('all')

class KappaR():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.rbins = np.geomspace(2e-4,32/128,32)
        self.r_cen= 0.5*(self.rbins[1:]+self.rbins[:-1])
    def run(self,core_list=None, do_plots=True, r_inflection=None, frame_list=None, tsing=None):
        this_looper=self.this_looper
        self.sim_name = this_looper.sim_name

        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        if frame_list is None:
            thtr=this_looper.tr
            mask = movie_frames.quantized_mask(this_looper).flatten()
            ok = np.zeros_like(mask)
            if 0:
                ok[::10] = mask[::10]
                mask=ok ;print('kludge mask')
            if 1:
                ok[-5]=True
                ok[-1]=True #just the last frame
                mask=ok
            times=thtr.times[mask]+0 #the zero makes a copy
            times.shape=times.size,1
            times=times/colors.tff
            self.times=times
            frame_list=thtr.frames[mask]
        self.frames = frame_list

        self.b_less = np.zeros([len(core_list),len(frame_list),len(self.rbins)])
        self.d_less = np.zeros([len(core_list),len(frame_list),len(self.rbins)])

        for ncore,core_id in enumerate(core_list):
            self.core_id=core_id
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            
            for nframe,frame in enumerate(frame_list):
                print("Flux on %s c%04d n%04d"%(this_looper.sim_name, core_id, frame))
                ds = this_looper.load(frame)
                nf = np.where( this_looper.tr.frames == frame)[0][0]
                time = times[nframe,0]

                c = nar([ms.mean_xc[nf], ms.mean_yc[nf],ms.mean_zc[nf]])
                msR = ms.rc
                msR[ msR<1/2048]=1/2048
                
                MaxRadius=msR[:,nf].max()
                Radius = max([8.0/128, MaxRadius])
                rsph = ds.arr(Radius,'code_length')
                sp = ds.sphere(c,rsph)

                R1 = sp['radius']
                order = np.argsort(R1)
                #vel = []
                #for axis in 'xyz':
                #    vel.append( sp['velocity_%s'%axis][order][:10].mean())
                #scrub = other_scrubber.scrubber(sp, reference_velocity = vel)
                #scrub.compute_ge()
                #scrub.compute_ke_rel()


                #Get data arrays

                #dv = scrub.cell_volume
                #RR = scrub.r
                #DD = scrub.density
                #divv = sp['velocity_divergence']
                dv = sp[YT_cell_volume]
                RR =sp[YT_radius]
                DD = sp[YT_density]
                Bi = sp[YT_magnetic_field_strength]

                ORDER = np.argsort( RR)
                RR_cuml = RR[ORDER]
                D_sort = DD[ORDER]
                dv_sort = dv[ORDER]
                M_local = D_sort*dv_sort
                M_cuml = np.cumsum( M_local )
                V_cuml = np.cumsum( dv[ORDER])
                D_cuml = M_cuml/V_cuml
                B_cuml = np.cumsum(Bi[ORDER]*M_local)/M_cuml


                digitized = np.digitize( RR_cuml, self.rbins)
                B_cuml_quant  = nar([ B_cuml[digitized==i].mean() if (digitized==i).any() else np.nan for i in range(0,len(self.rbins))])
                D_cuml_quant  = nar([ D_cuml[digitized==i].mean() if (digitized==i).any() else np.nan for i in range(0,len(self.rbins))])
                
                #ok = ~np.isnan(B_cuml_quant)
                self.b_less[ncore][nframe] = B_cuml_quant
                self.d_less[ncore][nframe] = D_cuml_quant



import track_loader as TL
sim_list=['u501','u502','u503']
sim_list=['u502']
TL.load_tracks(sim_list)

if 0:
    import tsing
    reload(tsing)
    if 'tsing_tool' not in dir():
        tsing_tool={}
        for ns,sim in enumerate(sim_list):
            obj=tsing.te_tc(TL.loops[sim])
            tsing_tool[sim]=obj
            tsing_tool[sim].run()
if 'things' not in dir():
    things={}

def kappa_r_time(thing):
    lnB = np.log10(thing.b_less)
    lnD = np.log10(thing.d_less)
    core_axis=0
    time_axis=1
    r_axis=2
    nancount_b=((~np.isnan(lnB)).sum(axis=core_axis))
    nancount_d=((~np.isnan(lnD)).sum(axis=core_axis))
    mean_lnB = np.nansum(lnB,axis=core_axis)/nancount_b
    mean_lnD = np.nansum(lnD,axis=core_axis)/nancount_d

    shapeB = list(mean_lnB.shape)
    shapeB.insert(core_axis,1)
    mean_lnB.shape = tuple(shapeB)
    mean_lnD.shape = tuple(shapeB)
    cross = np.nansum(((lnB-mean_lnB)*(lnD-mean_lnD)), axis=core_axis)
    var = np.nansum(((lnD-mean_lnD)**2),axis=core_axis)
    do_plots=True
    if 1:
        for nf,frame in enumerate(thing.frames):
            top = cross[nf,:]
            bot = var[nf,:]
            fig,axes=plt.subplots(1,1)
            ax0=axes
            ax0.plot( thing.rbins*colors.length_units_au,top/bot)
            ax0.set(xscale='log',xlabel='r [AU]', ylabel='kappa')
            fig.savefig('plots_to_sort/kappa_time_n%04d'%frame)
    if do_plots:
        for nf,frame in enumerate(thing.frames):
            for nr, r in enumerate(thing.rbins):
                fig,axes=plt.subplots(1,1)
                ax0=axes
                this_b = thing.b_less[:,nf,nr]
                this_d = thing.d_less[:,nf,nr]

                ax0.scatter(this_d,this_b)
                top = cross[nf,nr]
                bot=var[nf,nr]
                this_kappa=top/bot
                this_off  = mean_lnB[0,nf,nr]-this_kappa*mean_lnD[0,nf,nr]
                ln_d = np.log10(this_d)
                ax0.plot(this_d,10**(this_kappa*ln_d+this_off))
                ax0.set(xscale='log',yscale='log',xlabel='density',ylabel='total B')
                fig.savefig('plots_to_sort/many_plots_n%04d_r%04d'%(frame,nr))
            


        



if 1:
    #new plots.
    for sim in sim_list:
        if sim in things:
            continue
        all_cores=TL.loops[sim].core_list
        core_list = TL.loops[sim].core_by_mode['Alone']
        #core_list=core_list[4:8]

        #core_list = [114]
        thing = KappaR( TL.loops[sim])
        thing.run(core_list=core_list)
        things[sim]=thing

    if 1:
        for sim in sim_list:
            kappa_r_time(things[sim])
