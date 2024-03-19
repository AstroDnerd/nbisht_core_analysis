
# edit
from starter2 import *
import davetools
reload(davetools)
import p49_fields
reload(p49_fields)
import math
import pcolormesh_helper as pch
from matplotlib.ticker import PercentFormatter
np.set_printoptions(threshold=sys.maxsize)

import track_loader as TL
import cfpack as cfp
from cfpack import stop,print
import latex
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)

import tsing
reload(tsing)
import tsung_spheres
reload(tsung_spheres)
# --- --- --- --- --- --- ---


class withspheres(): 
    def __init__(self,the_loop):
        self.this_looper = the_loop
        self.cores_used = []
        self.bad_cores = []

        self.bmag_sph = defaultdict(list)
        self.rhoave_sph = defaultdict(list)
        self.bmag_parts = defaultdict(list)
        self.rhoave_parts = defaultdict(list)

    def framescores(self,sim,core_list=None): 
        print('inside!!')
        thtr = self.this_looper.tr

        # CORES
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores  #or debug

        # EVERY TEN FRAMES - make a loop now!
        for nf,frame in enumerate(thtr.frames): 

            # CORE-LOOP
            for nc,core_id in enumerate(core_list):
                self.cores_used.append(core_id)

                ds = self.this_looper.load(frame)  #loop frames 
                ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
                ms.particle_pos(core_id)

                # WITH SPHERES
                the_center = ms.mean_center[:,nf]  #the three coords for the last frame, will need to stop here. 
                the_radius = 1/128  #one zone   
                the_sphere = ds.sphere(the_center, the_radius) 
                # get and store the avg fields
                self.bmag_sph[nf].append((the_sphere['density']* the_sphere['magnetic_field_strength'] * the_sphere['cell_volume']).sum()/the_sphere['cell_mass'].sum()) 
                self.rhoave_sph[nf].append((the_sphere['density'] * the_sphere['cell_volume']).sum()/ the_sphere['cell_volume'].sum())

                # WITH PARTICLES 
                mask = ms.compute_unique_mask(core_id, dx=1./2048,frame=nf)
                # get the fields.
                density = thtr.c([core_id],'density')[mask,nf]
                cell_volume = thtr.c([core_id],'cell_volume')[mask,nf]
                cell_mass = density * cell_volume  
                bx = thtr.c([core_id],'magnetic_field_x')[mask,nf]
                by = thtr.c([core_id],'magnetic_field_y')[mask,nf]
                bz = thtr.c([core_id],'magnetic_field_z')[mask,nf]                
                bb = np.sqrt(bx*bx+by*by+bz*bz) 
                # store the avgs
                self.bmag_parts[nf].append((bb * density * cell_volume).sum()/cell_mass.sum())  
                self.rhoave_parts[nf].append((density * cell_volume).sum()/(cell_volume.sum()))  
        
        data_bmagsph = [*self.bmag_sph.values()] 
        data_rhosph = [*self.rhoave_sph.values()]
        data_bmagparts = [*self.bmag_parts.values()] 
        data_rhoparts = [*self.rhoave_parts.values()]
        if 0: 
            hfivename = 'p66_brho/brho_sphparts_%s.h5'%(sim)
            Fptr = h5py.File(hfivename,'w')
            Fptr['bfield_sph'] = data_bmagsph 
            Fptr['rhoavg_sph'] = data_rhosph
            Fptr['bfield_parts'] = data_bmagparts 
            Fptr['rhoavg_parst'] = data_rhoparts
            Fptr.close()
            print('h5 file written. closing.')



# YOU ENTER HERE
# TO GET DATA AND STORE
# COMPARING SPHERES WITH PARTICLES
if 0:
    sims=['u603', 'u602','u603']
    TL.load_tracks(sims)
    for sim in sims:
        core_list=None
        running = withspheres(TL.loops[sim])
        running.framescores(sim)
# SPHERES SYNCED TO TSUNG
if 0:
    sims=['u501', 'u502', 'u503']
    import three_loopers_u500 as TL   #EDIT THIS!! and put pdb.set_trace() back in this file
    #TL.load_tracks(sims)
    if 'tsing_tool' not in dir():
        tsing_tool={}
        for ns,sim in enumerate(sims):
            obj=tsing.te_tc(TL.loops[sim])
            tsing_tool[sim]=obj
            tsing_tool[sim].run()
    if 'mp' not in dir():
        for sim in sims:
            all_cores=np.unique(TL.loops[sim].tr.core_ids)
            core_list=list(all_cores)
            #core_list=core_list[:1]  #debug

            mp=tsung_spheres.tsungspheres(TL.loops[sim])   #EDIT! make tsung part of this class, then split THIS file into two respectively
            timescale = 2 
            mp.run(core_list=core_list, tsing=tsing_tool[sim], timescale=timescale, get_particles=True, save_sorts=True )



# TO READ DATA FROM STORAGE
if 1:  
    figtype = 'kappaperframe'  #kappaperframe, kappatff
    individually = 'no'  #yes, no
    parts_or_spheres ='sph_tsung' #parts, spheres, sphparts, sph_tsung  #EDIT respectively below

    sims=['u503']#, 'u502', 'u503']
    hfivename = 'p66_brho/h5files/brho_sphtsung_%s.h5'%(sims[0])  #EDIT
    Fptr = h5py.File(hfivename,'r')
    other_h5fields = ['bfield_sph','rhoavg_sph']
    b_sph = Fptr['bmag_sph'][()] 
    rho_sph = Fptr['rho_sph'][()]
    if parts_or_spheres == 'parts':
        b_parts = Fptr['bfield_parts'][()] 
        rho_parts = Fptr['rhoavg_parst'][()] 

    def afunct(x, a, b):
        y = a * x + b
        return y

    if figtype == 'kappaperframe':
        kappas_sph = []
        kappas_parts = []
        for i in range(len(rho_sph)):   #should make this part of the stored data... 
            rhosph_log = np.log(rho_sph[i])  
            bsph_log = np.log(b_sph[i]) 
            rets_sph = cfp.fit(afunct, rhosph_log, bsph_log)  
            kappas_sph = np.append(kappas_sph, rets_sph.popt[0])
            if parts_or_spheres == 'parts':
                bparts_log = np.log(b_parts[i]) 
                rhoparts_log = np.log(rho_parts[i])  
                rets_parts = cfp.fit(afunct, rhoparts_log, bparts_log)  
                kappas_parts = np.append(kappas_parts, rets_parts.popt[0]) 
            '''
            the_xrecipe = np.linspace(rhosph_log.min(),rhosph_log.max(),num=500)
            the_xten =10**the_xrecipe
            the_yrecipe_ = afunct(the_xrecipe, *rets.popt)
            the_yten = 10**the_yrecipe
            ''' 
            if individually == 'yes':
                fig,ax = plt.subplots(1,1)
                ax.scatter(rho_sph[i], b_sph[i], color='b', label='spheres', alpha=0.5)  #kappa sph and parts vs time frame
                if parts_or_spheres == 'parts':
                    ax.scatter(rho_parts[i], b_parts[i], color='r', label='particles', alpha=0.5)  #kappa sph and parts vs time frame
                ax.legend()
                ax.set(xlabel=r'$\rho_{ave}$', ylabel=r'$|B|$', xscale='log', yscale='log', xlim=(1e-2,1e8), ylim=(1e0,1e4))
                outname = 'p66_brho/sphtsung_frame%d_%s'%(i,sims[0])
                plt.savefig(outname)
                print('figure saved!')
                plt.clf()      
        if individually == 'no':
            fig,ax = plt.subplots(1,1)
            the_x = np.linspace(1,len(rho_sph),len(rho_sph)) 
            ax.scatter(the_x, kappas_sph, color='b', label='spheres')  #kappa sph and parts vs time frame
            if parts_or_spheres == 'parts':
                ax.plot(the_x, kappas_parts, color='r', label='particles')  #kappa sph and parts vs time frame
            ax.legend()
            ax.set(xlabel=r'$t_{tsung,dummy}$', ylabel=r'$\kappa_{sph_tsung}$') #,\kappa_{parts}$', ylim=(0,0.95))
            outname = 'p66_brho/sphtsung_kappadyn_scatter_%s'%sims[0]
            plt.savefig(outname)
            print('figure saved!')
            plt.clf()    

    if figtype == 'kappatff':
        rhotff_sph = []
        btff_sph = []
        rhotff_parts = []
        btff_parts = []
        fig,ax = plt.subplots(1,1)
        for i in range(len(rho_sph)):   
            rhotff_sph = np.append(rhotff_sph, rho_sph[i]) 
            btff_sph = np.append(btff_sph, b_sph[i]) 
            if parts_or_spheres == 'parts':
                rhotff_parts = np.append(rhotff_parts, rho_parts[i]) 
                btff_parts = np.append(btff_parts, b_parts[i])   

        rhotff_sphlog = np.log10(rhotff_sph)
        btff_sphlog = np.log10(btff_sph) 
        rets_sphlog = cfp.fit(afunct, rhotff_sphlog, btff_sphlog)  
        if parts_or_spheres == 'parts':
            rhotff_partslog = np.log10(rhotff_parts)
            btff_partslog = np.log10(btff_parts)
            rets_partslog = cfp.fit(afunct, rhotff_partslog, btff_partslog)  

        tmap = rainbow_map(len(rhotff_sph)) 
        ctr = [tmap(n) for n in range(len(rhotff_sph))]
        color_opt = [ctr,'b','r']       
        if 1:
            ax.scatter(rhotff_sph, btff_sph, color=color_opt[0], alpha=0.4, label='spheres')
        if 0:
            ax.scatter(rhotff_parts, btff_parts, color=color_opt[0], alpha=0.4, label='particles')

        title=[r'$\kappa_{sph_tsung}=%f$'%rets_sphlog.popt[0]]#, r'$\kappa_{sph}=%f$, $\kappa_{parts}=%f$'%(rets_sphlog.popt[0],rets_partslog.popt[0]), 'particles', 'spheres']
        ax.set(xlabel=r'$\rho_{ave}$', ylabel=r'$|B|$', xscale='log', yscale='log', xlim=(10e-3,10e3), ylim=(10e-2,10e3), \
               title=title[0])  
        ax.legend(loc='best')

        outname = 'p66_brho/%s_kappatff_scatter_%s'%(parts_or_spheres,sims[0])
        plt.savefig(outname)
        print('figure saved!')
        plt.clf()    

    print('closing fig and h5 file!')
    Fptr.close()

# simple rho was a definition, withspheres is a class
# UFFF
