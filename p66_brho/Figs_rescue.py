
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
        if 1: 
            hfivename = 'p66_brho/brho_sphparts_%s.h5'%(sim)
            Fptr = h5py.File(hfivename,'w')
            Fptr['bfield_sph'] = data_bmagsph 
            Fptr['rhoavg_sph'] = data_rhosph
            Fptr['bfield_parts'] = data_bmagparts 
            Fptr['rhoavg_parst'] = data_rhoparts
            Fptr.close()
            print('h5 file written. closing.')


# YOU ARE HERE
sims=['u603']#, 'u602','u603']
# TO GET DATA AND STORE
if 0:
    TL.load_tracks(sims)
    for sim in sims:
        core_list=None
        running = withspheres(TL.loops[sim])
        running.framescores(sim)

# TO READ DATA FROM STORAGE
if 1:  
    figtype = 'kappaperframe'  #kappaperframe, kappatff
    individually = 'yes'  #yes, no
    parts_or_spheres ='parts' #parts, spheres, sphparts  #EDIT respectively below

    hfivename = 'p66_brho/brho_sphparts_%s.h5'%(sims[0])  #EDIT
    Fptr = h5py.File(hfivename,'r')
    b_sph = Fptr['bfield_sph'][()] 
    b_parts = Fptr['bfield_parts'][()] 
    rho_sph = Fptr['rhoavg_sph'][()]
    rho_parts = Fptr['rhoavg_parst'][()] 
    #pdb.set_trace()

    def afunct(x, a, b):
        y = a * x + b
        return y

    if figtype == 'kappaperframe':
        kappas_sph = []
        kappas_parts = []
        for i in range(len(rho_parts)):   #should make this part of the stored data... 
            rhosph_log = np.log(rho_sph[i])  
            rhoparts_log = np.log(rho_parts[i])  
            bsph_log = np.log(b_sph[i]) 
            bparts_log = np.log(b_parts[i]) 

            rets_sph = cfp.fit(afunct, rhosph_log, bsph_log)  
            rets_parts = cfp.fit(afunct, rhoparts_log, bparts_log)  
            '''
            the_xrecipe = np.linspace(rhosph_log.min(),rhosph_log.max(),num=500)
            the_xten =10**the_xrecipe
            the_yrecipe_ = afunct(the_xrecipe, *rets.popt)
            the_yten = 10**the_yrecipe
            '''
            kappas_sph = np.append(kappas_sph, rets_sph.popt[0])
            kappas_parts = np.append(kappas_parts, rets_parts.popt[0]) 
            
            if individually == 'yes':
                fig,ax = plt.subplots(1,1)
                ax.scatter(rho_sph[i], b_sph[i], color='b', label='spheres', alpha=0.5)  #kappa sph and parts vs time frame
                ax.scatter(rho_parts[i], b_parts[i], color='r', label='particles', alpha=0.5)  #kappa sph and parts vs time frame
                ax.legend()
                ax.set(xlabel=r'$\rho_{ave}$', ylabel=r'$|B|$', xscale='log', yscale='log', xlim=(1e-2,1e8), ylim=(1e0,1e4))
                outname = 'p66_brho/sphparts_frame%d_%s'%(i,sims[0])
                plt.savefig(outname)
                print('figure saved!')
                plt.clf()      
        if individually == 'no':
            fig,ax = plt.subplots(1,1)
            the_x = np.linspace(0,len(rho_parts)+1,len(rho_parts)) 
            ax.plot(the_x, kappas_sph, color='b', label='spheres')  #kappa sph and parts vs time frame
            ax.plot(the_x, kappas_parts, color='r', label='particles')  #kappa sph and parts vs time frame
            ax.legend()
            ax.set(xlabel=r'$t_{ff,dummy}$', ylabel=r'$\kappa_{sph},\kappa_{parts}$', ylim=(0,0.95))
            outname = 'p66_brho/sphparts_kappadyn_plot_%s'%sims[0]
            plt.savefig(outname)
            print('figure saved!')
            plt.clf()    

    if figtype == 'kappatff':
        rhotff_sph = []
        btff_sph = []
        rhotff_parts = []
        btff_parts = []
        fig,ax = plt.subplots(1,1)
        for i in range(len(rho_parts)):   
            rhotff_sph = np.append(rhotff_sph, rho_sph[i]) 
            btff_sph = np.append(btff_sph, b_sph[i]) 
            rhotff_parts = np.append(rhotff_parts, rho_parts[i]) 
            btff_parts = np.append(btff_parts, b_parts[i])   

        rhotff_sphlog = np.log10(rhotff_sph)
        btff_sphlog = np.log10(btff_sph) 
        rhotff_partslog = np.log10(rhotff_parts)
        btff_partslog = np.log10(btff_parts)
        rets_sphlog = cfp.fit(afunct, rhotff_sphlog, btff_sphlog)  
        rets_partslog = cfp.fit(afunct, rhotff_partslog, btff_partslog)  

        tmap = rainbow_map(len(rhotff_parts)) 
        ctr = [tmap(n) for n in range(len(rhotff_parts))]
        color_opt = [ctr,'b','r']
       
        if 0:
            ax.scatter(rhotff_sph, btff_sph, color=color_opt[0], alpha=0.4, label='spheres')
        if 1:
            ax.scatter(rhotff_parts, btff_parts, color=color_opt[0], alpha=0.4, label='particles')

        title=[r'$\kappa_{sph}=%f$, $\kappa_{parts}=%f$'%(rets_sphlog.popt[0],rets_partslog.popt[0]), 'particles', 'spheres']
        ax.set(xlabel=r'$\rho_{ave}$', ylabel=r'$|B|$', xscale='log', yscale='log', \
               title=title[1])  
        ax.legend(loc='best')

        outname = 'p66_brho/%s_kappatff_scatter_%s'%(parts_or_spheres,sims[0])
        plt.savefig(outname)
        print('figure saved!')
        plt.clf()    

    print('closing fig and h5 file!')
    Fptr.close()

# simple rho was a definition, withspheres is a class
# UFFF
