
'''
needs to be in the same folder as tsung_spheres
'''
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
import lmfit
import latex
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import seaborn as sns

import tsing
reload(tsing)
import tsung_spheres
reload(tsung_spheres)
import monster
reload(monster)

from yt.data_objects.level_sets.api import *
# --- --- --- --- --- --- ---

class withspheres(): 
    def __init__(self,the_loop):
        self.this_looper = the_loop

        self.bmag_sph_rsmtwo = defaultdict(list)
        self.rhoave_sph_rsmtwo = defaultdict(list)

        self.bmag_smparts = defaultdict(list)
        self.rhoave_smparts = defaultdict(list)

        self.ncolumn_sph_rsmtwo = defaultdict(list)
        self.blos_sph_rsmtwo = defaultdict(list)
        

    def framescores(self, sim, typeofstudy='withspheres', core_list=None): 
        print('inside!!')
        monster.load([sim])
        the_monster = monster.closet[sim]
        thtr = the_monster.this_looper.tr  #NEW - check,good

        # CORES
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            print('all cores')
            corelist = all_cores
            #corelist = all_cores[1:2]  #or debug
        if core_list == 'alone':
            print('lonely cores')
            corelist = the_monster.this_looper.core_by_mode['A']
            #corelist=corelist[1:2]

        # EVERY TEN FRAMES
        for nf,frame in enumerate(thtr.frames[1:]):  #or debug 

            # CORE-LOOP
            for nc,core_id in enumerate(corelist):
                the_ms = the_monster.get_ms(core_id, do_velocity=False, do_magnetic=True)  #NEW

                # get and store the avg fields
                # WITH SPHERES, R_SMART_2, WEIGHT=MASS
                the_sphere_rsmtwo = the_monster.get_sphere(core_id,frame,'rsmart_2')
                mass_total = the_sphere_rsmtwo['cell_mass'].sum()
                volume_total = the_sphere_rsmtwo['cell_volume'].sum()
                self.bmag_sph_rsmtwo[nf].append((the_sphere_rsmtwo['cell_mass']* the_sphere_rsmtwo['magnetic_field_strength']).sum()/mass_total) 
                self.rhoave_sph_rsmtwo[nf].append(mass_total/volume_total)

                # WITH PARTICLES 
                if typeofstudy == 'particles':
                    mask = the_ms.compute_unique_mask(core_id, dx=1./2048,frame=nf)
                    # get the fields.
                    density = thtr.c([core_id],'density')[mask,nf]
                    cell_volume = thtr.c([core_id],'cell_volume')[mask,nf]
                    cell_mass = density*cell_volume
                    bx = thtr.c([core_id],'magnetic_field_x')[mask,nf]
                    by = thtr.c([core_id],'magnetic_field_y')[mask,nf]
                    bz = thtr.c([core_id],'magnetic_field_z')[mask,nf]                
                    bb = np.sqrt(bx*bx+by*by+bz*bz) 
                    # store the avgs
                    self.bmag_smparts[nf].append((bb * cell_mass).sum()/cell_mass.sum())  
                    self.rhoave_smparts[nf].append(cell_mass.sum()/cell_volume.sum())  
     
        if typeofstudy == 'particles':
            data_bmagparts = [*self.bmag_smparts.values()] 
            data_rhoparts = [*self.rhoave_smparts.values()]
        data_bmagsph = [*self.bmag_sph_rsmtwo.values()] 
        data_rhosph = [*self.rhoave_sph_rsmtwo.values()]

        print('before saving to h5 files!!')
        if 1: 
            hfivename = 'p66_brho/brho_sphparts_rsmtwo_all_%s.h5'%(sim)
            Fptr = h5py.File(hfivename,'w')
            Fptr['bfield_sph'] = data_bmagsph 
            Fptr['rhoavg_sph'] = data_rhosph
            if typeofstudy == 'particles':
                Fptr['bfield_parts'] = data_bmagparts 
                Fptr['rhoavg_parts'] = data_rhoparts
            Fptr.close()
            print('h5 file written. closing.')


    def syntheticobs(self, sim, core_list=None): 
        print('inside synthetic observations!')
        monster.load([sim])
        the_monster = monster.closet[sim]
        thtr = the_monster.this_looper.tr  

        # CORES
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            corelist = all_cores
            #corelist = all_cores[1:2] #or debug
        if core_list == 'alone':
            corelist = the_monster.this_looper.core_by_mode['A']

        # EVERY TEN FRAMES 
        for nf,frame in enumerate(thtr.frames[1:]): 
            # CORE-LOOP
            for nc,core_id in enumerate(corelist):
                the_sphere_rsmtwo = the_monster.get_sphere(core_id,frame,'rsmart_2') 
                the_rsmtwo_rad = the_monster.get_rsmarttwo(core_id,frame) 
                the_area_rsmtwo = np.pi * the_rsmtwo_rad**2
                
                B = ['magnetic_field_x','magnetic_field_y','magnetic_field_z']
                for j in range(3): 
                    self.ncolumn_sph_rsmtwo[nf].append(the_sphere_rsmtwo['cell_mass'].sum()/the_area_rsmtwo)
                    self.blos_sph_rsmtwo[nf].append((the_sphere_rsmtwo['cell_mass'] * the_sphere_rsmtwo[B[j]]).sum()/the_sphere_rsmtwo['gas','cell_mass'].sum())
  
        data_blossph = [*self.blos_sph_rsmtwo.values()] 
        data_ncolumnsph = [*self.ncolumn_sph_rsmtwo.values()]
        print('before saving to h5 files')
        if 1: 
            hfivename = 'p66_brho/blosncol_rsmtwo_alone_%s.h5'%(sim)
            Fptr = h5py.File(hfivename,'w')
            Fptr['blos_sph'] = data_blossph 
            Fptr['ncolumn_sph'] = data_ncolumnsph
            Fptr.close()
            print('h5 file written. closing.')



    def projspheres(self, sim, core_list=None, individual='no'): 
        print('inside projections')
        #thtr = self.this_looper.tr
        monster.load([sim])
        the_monster = monster.closet[sim]
        thtr = the_monster.this_looper.tr  #NEW - check

        # CORES
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            print('all cores')
            corelist = all_cores #or debug
        if core_list == 'alone':
            print('lonely cores')
            corelist = the_monster.this_looper.core_by_mode['A']

        # EVERY TEN FRAMES
        for nf,frame in enumerate(thtr.frames[1:]):  #or debug 
            ds = the_monster.get_ds(frame) 
            if individual == 'no':
                p = yt.ProjectionPlot(ds, 2, ("gas", "density"))

            # CORE-LOOP
            for nc,core_id in enumerate(corelist):
                #if core_id == 109:
                if 1: 
                    the_sphere_rinf = the_monster.get_sphere(core_id,frame,'rinf')
                    the_sphere_rone = the_monster.get_sphere(core_id,frame,'r1')
                    the_sphere_reight = the_monster.get_sphere(core_id,frame,'r8')
                    the_sphere_rmax = the_monster.get_sphere(core_id,frame,'rmax')

                    the_radius_rinf = the_monster.get_r_inflection(core_id,frame)
                    the_radius_eight = 8/128
                    the_radius_one = 1/128
                    the_center = the_sphere_rinf.center 

                    if individual == 'no':
                        p.annotate_sphere(the_center, radius=the_radius_rinf, circle_args={"color": "red"}, text='%d'%core_id)
                        if 0:
                            p.annotate_sphere(the_center, radius=the_radius_one, circle_args={"color": "black"}, text='%d'%core_id)
                    if individual == 'yes':
                        p = yt.ProjectionPlot(ds, 2, ("gas", "density"), center=the_center, data_source=the_sphere_reight)
                        p.set_width(2*the_radius_eight)
                        p.annotate_sphere(the_center, radius=the_radius_rinf, circle_args={"color": "red"})
                        p.annotate_sphere(the_center, radius=the_radius_one, circle_args={"color": "black"})
                        print('about to save')
                        p.save('./p66_brho/%d_%d_%s'%(core_id,frame,sim))

            if individual == 'no':
                print('saving proj')
                p.save('./p66_brho/%d_%s'%(frame,sim))



# ENTER HERE
# TO GET DATA AND STORE
# COMPARING SPHERES WITH PARTICLES
# note: take a look at p19_play/pseudo.py for a brief example with monster.py
if 0:
    sims=['u601']#, 'u602','u603']
    TL.load_tracks(sims)
    for sim in sims:
        the_corelist='alone'
        running = withspheres(TL.loops[sim])
        if 0:
            running.framescores(sim,typeofstudy='particles',core_list=the_corelist)
        if 1:
            running.syntheticobs(sim,core_list=the_corelist)
        if 0:
            running.projspheres(sim, individual='no')
# TSUNG TIMES
if 0:
    sims=['u501', 'u502', 'u503']
    TL.load_tracks(sims)
    for sim in sims:
        monster.load([sim])
        the_monster = monster.closet[sim]
        thtr = the_monster.this_looper.tr 
        all_cores = np.unique(thtr.core_ids)
        if 1:
            core_list = all_cores  #or debug
        if 0:
            core_list = the_monster.this_looper.core_by_mode['A']
        tsungtff =[]
        alltime = thtr.times
        for nc,core_id in enumerate(core_list):
            tsung_time = the_monster.get_tsung(core_id) 
            tsungtff.append(thtr.times[the_monster.get_time_index(tsung_time)])
        stop()

G = 1620/(4*np.pi)
rho_mean = 1
t_ff = np.sqrt(3*np.pi/(32*G*rho_mean))




# TO READ DATA FROM STORAGE
if 1:   
    whichcores = 'alone'
    spheres = 'yes'
    synthetic = 'no'
    particles = 'no'
    figtype = 'kappatff_framed'  #kappaperframe, kappatff, kappatff_framed
    individually = 'no'  #if kappaperframe, yes: one panel per frame, no: frame time series
    corestory = 'no'    #if kappatff, yes: one core per tff, no: all cores per tff
    save='no'
    tff_core='no'

    if spheres =='yes' and particles=='yes':
        saving='sphparts'
    if spheres =='yes' and particles=='no':
        saving='sph'
    if synthetic =='yes' and particles=='no':
        saving='synth'
    if spheres =='no' and particles=='yes':
        saving='parts'
    if synthetic =='yes' and spheres=='yes':
        saving='sphsynth'
    if spheres =='yes' and corestory=='yes':
        saving='sphsing'

    sims=['u603']#, 'u602', 'u603']  # need a loop!!
    hfivename = 'p66_brho/h5files/brho_sphparts_rsmtwo_alone_%s.h5'%(sims[0])  
    Fptr = h5py.File(hfivename,'r')
    b_sph = Fptr['bfield_sph'][()] 
    rho_sph = Fptr['rhoavg_sph'][()] 
    b_parts = Fptr['bfield_parts'][()] 
    rho_parts = Fptr['rhoavg_parts'][()] 
    if synthetic == 'yes':
        hfivename_synth = 'p66_brho/h5files/blosncol_rsmtwo_alone_%s.h5'%(sims[0])  
        Fptr_synth = h5py.File(hfivename_synth,'r')
        b_sph_synth = Fptr_synth['blos_sph'][()] 
        rho_sph_synth = Fptr_synth['ncolumn_sph'][()] 

    def afunct(x, a, b):
        y = a * x + b
        return y


    if figtype == 'kappaperframe':
        kappas_sph = []
        kappas_parts = []
        kappas_sph_synth = []
        for i in range(len(rho_sph)):   
            if spheres == 'yes':
                rhosph_log = np.log10(rho_sph[i])
                bsph_log = np.log10(abs(b_sph[i]))
                rets_sph = cfp.fit(afunct, rhosph_log, bsph_log)  
                kappas_sph = np.append(kappas_sph, rets_sph.popt[0])
            if synthetic == 'yes':
                rhosph_synth_log = np.log10(rho_sph_synth[i])
                bsph_synth_log = np.log10(abs(b_sph_synth[i]))
                rets_sph_synth = cfp.fit(afunct, rhosph_synth_log, bsph_synth_log)  
                kappas_sph_synth = np.append(kappas_sph_synth, rets_sph_synth.popt[0])
            if particles == 'yes':
                bparts_log = np.log10(b_parts[i]) 
                rhoparts_log = np.log10(rho_parts[i])  
                rets_parts = cfp.fit(afunct, rhoparts_log, bparts_log)  
                kappas_parts = np.append(kappas_parts, rets_parts.popt[0]) 

            # one panel for each time frame
            if individually == 'yes':
                fig,ax = plt.subplots(1,1)
                if spheres == 'yes':
                    ax.scatter(rho_sph[i], b_sph[i], color='b', label='rsmart2', alpha=0.5)  
                    the_xrecipe = np.linspace(rhosph_log.min(),rhosph_log.max(),num=500)
                    the_yrecipe = afunct(the_xrecipe, *rets_sph.popt)
                    the_xten =10**the_xrecipe
                    the_yten = 10**the_yrecipe
                    ax.plot(the_xten, the_yten, c='gray', linewidth=2.0, alpha=0.7)
                if synthetic == 'yes':
                    ax.scatter(rho_sph_synth[i], b_sph_synth[i], color='r', label='synthetic', alpha=0.5)  
                    the_xrecipe = np.linspace(rhosph_synth_log.min(),rhosph_synth_log.max(),num=500)
                    the_yrecipe = afunct(the_xrecipe, *rets_sph_synth.popt)
                    the_xten_synth =10**the_xrecipe
                    the_yten_synth = 10**the_yrecipe
                    ax.plot(the_xten_synth, the_yten_synth, c='gray', linewidth=2.0, alpha=0.7)
                if particles == 'yes':
                    ax.scatter(rho_parts[i], b_parts[i], color='r', label='particles', alpha=0.5)  
                    the_xrecipe = np.linspace(rparts_log.min(),rhoparts_log.max(),num=500)
                    the_yrecipe = afunct(the_xrecipe, *rets_parts.popt)
                    the_xten_parts =10**the_xrecipe
                    the_yten_parts = 10**the_yrecipe
                    ax.plot(the_xten_parts, the_yten_parts, c='gray', linewidth=2.0, alpha=0.7)

                ax.legend(loc='best')
                #title=r'$\kappa_{sph}=%f, \kappa_{parts}=%f $'%(rets_sph.popt[0],rets_parts.popt[0])
                ax.set(xlabel=r'$\rho_{ave}$', ylabel=r'$|B|$', xscale='log', yscale='log', xlim=(1e-2,1e8), ylim=(1.5e-2,1.5e4), \
                    title=None)
                outname = 'p66_brho/rsmartparts_frame%d_%s'%(i,sims[0])
                plt.savefig(outname)
                print('figure saved!')
                plt.clf()      
                
        # time frames in time series
        if individually == 'no':
            fig,ax = plt.subplots(1,1)
            the_x = np.linspace(1,len(rho_sph),len(rho_sph))  
            if spheres == 'yes':
                ax.scatter(the_x, kappas_sph, color='k', label='spheres')  
            if synthetic == 'yes':
                ax.scatter(the_x, kappas_sph_synth, color='g', label='synthetic')  
            if particles == 'yes':
                ax.scatter(the_x, kappas_parts, color='r', label='particles')  

            outname = 'p66_brho/%s_kappadyn_scatter_%s_%s'%(saving,whichcores,sims[0]) 
            ax.legend(loc='best')
            ax.set(ylabel=r'$\kappa$', xlabel=r'$t_{tff,dummy}$', ylim=(0,1.0))
            plt.savefig(outname)
            print('figure saved!')
            plt.clf()    


    if figtype == 'kappatff':
        rhotff_sph = []
        btff_sph = []
        rhotff_parts = []
        btff_parts = []
        rhotff_synth = []
        btff_synth = []
        fig,ax = plt.subplots(1,1)

        if corestory=='no': 
            for i in range(len(rho_sph)):   
                if spheres == 'yes':
                    rhotff_sph = np.append(rhotff_sph, rho_sph[i]) 
                    btff_sph = np.append(btff_sph, b_sph[i]) 
                if particles == 'yes':
                    rhotff_parts = np.append(rhotff_parts, rho_parts[i]) 
                    btff_parts = np.append(btff_parts, b_parts[i])   
                if synthetic == 'yes':
                    rhotff_synth = np.append(rhotff_synth, rho_sph_synth[i]) 
                    btff_synth = np.append(btff_synth, b_sph_synth[i])   
            if spheres == 'yes':
                rhotff_sphlog = np.log10(rhotff_sph)
                btff_sphlog = np.log10(btff_sph) 
                rets_sphlog = cfp.fit(afunct, rhotff_sphlog, btff_sphlog)  
                tmap = rainbow_map(len(rhotff_sph)) 
                ctr = [tmap(n) for n in range(len(rhotff_sph))]
            if synthetic == 'yes':
                rhotff_synthrinflog = np.log10(rhotff_synth)
                btff_synthrinflog = np.log10(abs(btff_synth))
                rets_synthlog = cfp.fit(afunct, rhotff_synthrinflog, btff_synthrinflog)  
                tmap = rainbow_map(len(rhotff_synth)) 
                ctr = [tmap(n) for n in range(len(rhotff_synth))]
            if particles == 'yes':
                rhotff_partslog = np.log10(rhotff_parts)
                btff_partslog = np.log10(btff_parts)
                rets_partslog = cfp.fit(afunct, rhotff_partslog, btff_partslog)  
                tmap = rainbow_map(len(rhotff_parts)) 
                ctr = [tmap(n) for n in range(len(rhotff_parts))]

            color_opt = [ctr,'b','r']       
            if spheres == 'yes':
                ax.scatter(rhotff_sph, btff_sph, color=color_opt[0], alpha=0.4, label='spheres')
                the_xrecipe = np.linspace(rhotff_sphlog.min(),rhotff_sphlog.max(),num=500)
                the_yrecipe = afunct(the_xrecipe, *rets_sphlog.popt)
                the_xten =10**the_xrecipe
                the_yten = 10**the_yrecipe
                ax.plot(the_xten, the_yten, c='gray', linewidth=2.0, alpha=0.7)
                title=r'$\kappa_{sph}=:%f$'%rets_sphlog.popt[0]
            if synthetic == 'yes':
                ax.scatter(rhotff_synth, btff_synth, color=color_opt[0], alpha=0.4, label='synth')
                the_xrecipe = np.linspace(rhotff_synthrinflog.min(),rhotff_synthrinflog.max(),num=500)
                the_yrecipe = afunct(the_xrecipe, *rets_synthlog.popt)
                the_xten =10**the_xrecipe
                the_yten = 10**the_yrecipe
                ax.plot(the_xten, the_yten, c='gray', linewidth=2.0, alpha=0.7)
                title=r'$\kappa_{synth}=%f$'%rets_synthlog.popt[0]
            if particles == 'yes':
                if spheres =='yes':
                    ax.scatter(rhotff_parts, btff_parts, color='k', alpha=0.4, label='particles')
                else:
                    ax.scatter(rhotff_parts, btff_parts, color=color_opt[0], alpha=0.4, label='particles')
                the_xrecipe = np.linspace(rhotff_partslog.min(),rhotff_partslog.max(),num=500)
                the_yrecipe = afunct(the_xrecipe, *rets_partslog.popt)
                the_xten =10**the_xrecipe
                the_yten = 10**the_yrecipe
                ax.plot(the_xten, the_yten, c='gray', linewidth=2.0, alpha=0.7)
                title=r'$\kappa_{parts}=%f$'%rets_partslog.popt[0]

            if particles =='yes' and spheres=='yes':
                title = None

            outname = 'p66_brho/%s_kappatff_scatter_%s_%s'%(saving, whichcores, sims[0])
            ax.set(xlabel=r'$\rho_{ave}$', ylabel=r'$|B|$', xscale='log', yscale='log', title=title)
            ax.legend(loc='best')
            plt.savefig(outname)
            print('figure saved!')
            plt.clf()    


        if corestory=='yes':
            kappas=[]
            for j in range(len(rho_sph[0])): 
                for i in range(len(rho_sph)):
                    if spheres == 'yes':
                        rhotff_sph = np.append(rhotff_sph, rho_sph[i][j]) 
                        btff_sph = np.append(btff_sph, b_sph[i][j]) 
                if spheres == 'yes':
                    rhotff_sphlog = np.log10(rhotff_sph)
                    btff_sphlog = np.log10(btff_sph) 
                    rets_sphlog = cfp.fit(afunct, rhotff_sphlog, btff_sphlog)  
                    tmap = rainbow_map(len(rhotff_sph)) 
                    ctr = [tmap(n) for n in range(len(rhotff_sph))]
                color_opt = [ctr,'b','r']       
                if spheres == 'yes':
                    the_xrecipe = np.linspace(rhotff_sphlog.min(),rhotff_sphlog.max(),num=500)
                    the_yrecipe = afunct(the_xrecipe, *rets_sphlog.popt)
                    the_xten =10**the_xrecipe
                    the_yten = 10**the_yrecipe
                    title=r'$\kappa_{sph}=:%f$'%rets_sphlog.popt[0]
                kappas.append(rets_sphlog.popt[0])
                if save=='yes':
                    ax.scatter(rhotff_sph, btff_sph, color=color_opt[0], alpha=0.4, label='spheres')
                    ax.plot(the_xten, the_yten, c='gray', linewidth=2.0, linestyle='dotted', alpha=0.5)
                    ax.set(xlabel=r'$\rho_{ave}$', ylabel=r'$|B|$', xscale='log', yscale='log', title=title)
                    ax.legend(loc='best')
                    outname = 'p66_brho/%s_kappatff_scatter_%s_%s'%(saving, whichcores, sims[0])
                    plt.savefig(outname)
                    print('figure saved!')
                    plt.clf()    

            # histogram of tff per core!
            if tff_core=='yes':
                kappas = np.array(kappas)
                the_min = kappas.min()
                the_max = kappas.max()
                the_bins = np.linspace(the_min, the_max, num=64)  #there was a comment made at the group meet about this.
                ax.hist(kappas, bins=the_bins, density=True, histtype='step', color='k')
                kappas_mean = kappas.mean()
                ax.set(title='mean %s core per tff kappa: %f'%(whichcores, kappas_mean), xlim=(0.4,0.7), \
                       xlabel=r'$\kappa_{\mathrm{core_tff}}$', ylabel='count')
                print('saving histos')
                plt.savefig('p66_brho/kappastff_percore_%s_%s'%(whichcores, sims[0]))



    if figtype=='kappatff_framed':                
        fig,ax = plt.subplots(1,1)
        kappas_sph=[]
        rhotff_sph=[]
        btff_sph=[]
        for j in range(len(rho_sph)):   
            rhotff_sph = np.append(rhotff_sph, rho_sph[j]) 
            btff_sph = np.append(btff_sph, b_sph[j]) 
        rhotff_sphlog = np.log10(rhotff_sph)
        btff_sphlog = np.log10(btff_sph) 
        rets_sphlog = cfp.fit(afunct, rhotff_sphlog, btff_sphlog)  

        tmap = rainbow_map(len(rhotff_sph)) 
        ctr = [tmap(n) for n in range(len(rhotff_sph))]
        color_opt = [ctr,'b','r']       
        ax.scatter(rhotff_sph, btff_sph, color=color_opt[0], alpha=0.4)#, label='spheres')

        for i in range(len(rho_sph)):   
            rhosph_log = np.log10(rho_sph[i])
            bsph_log = np.log10(abs(b_sph[i]))
            rets_sph = cfp.fit(afunct, rhosph_log, bsph_log)  
            kappas_sph = np.append(kappas_sph, rets_sph.popt[0])
            the_xrecipe = np.linspace(rhosph_log.min(),rhosph_log.max(),num=500)
            the_yrecipe = afunct(the_xrecipe, *rets_sph.popt)
            the_xten =10**the_xrecipe
            the_yten = 10**the_yrecipe
            ax.plot(the_xten, the_yten, color='k', linewidth=2.0, alpha=0.4)

        the_xrecipe = np.linspace(rhotff_sphlog.min(),rhotff_sphlog.max(),num=500)
        the_yrecipe = afunct(the_xrecipe, *rets_sphlog.popt)
        the_xten =10**the_xrecipe
        the_yten = 10**the_yrecipe

        #ax.plot(the_xten, the_yten, c='gray', linewidth=2.0, alpha=0.7)
        #title=r'$\kappa_{sph}=:%f$'%rets_sphlog.popt[0]
        title=None
        outname = 'p66_brho/%s_kappatffframes_scatter_%s_%s'%(saving, whichcores, sims[0])
        ax.set(xlabel=r'$\rho_{ave}$', ylabel=r'$|B|$', xscale='log', yscale='log', title=title)
        #ax.legend(loc='best')
        plt.savefig(outname)
        print('figure saved!')
        plt.clf()    
        


    # tsun of t_ff, keep for now
    if figtype == 'tsungtff':
        fig,ax = plt.subplots(1,1)
        tff_ratio = tsung_tff/t_ff
        the_min = tff_ratio.min()
        the_max = tff_ratio.max()
        the_bins = np.linspace(the_min, the_max, num=64)  #there was a comment made at the group meet about this.
        ax.hist(tff_ratio, bins=the_bins, density=True, histtype='step', color='k')
        tsungtffratio = tff_ratio.mean()

        ax.set(title='mean tsung/tff ratio: %f'%tsungtffratio, xlabel=r'$t_{\mathrm{sung}}/t_{\mathrm{ff}}$', ylabel='count')
        print('saving histos')
        plt.savefig('p66_brho/plotexplore/tsungratiotff_%s'%sims[0])

    print('closing fig and h5 file!')
    Fptr.close()





