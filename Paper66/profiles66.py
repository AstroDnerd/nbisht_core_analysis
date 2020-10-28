
"""
 PART OF P19_NEWSCRIPTS REPOSITORY:
This script must be placed on: ~/p19_newscripts
With: starter1.py, starter2.py on same directory.
Started as: ~/p19_newscripts/examples/mask_profile.py
In conjuction with: Brho_particles66.py


A SCRIPT FOR PROFILES,PHASES,PLOT
NOTE on DATA:
/scratch2/luzlourdes/simulations/u05grav/GravPotential, faulty disk, careful
"""
from starter2 import *
from starter1 import *
import Brho_particles66 as brp 


def _BoverRho(field,data):
    return data['magnetic_field_strength']/data['density']
yt.add_field('_BoverRho',function=_BoverRho,units='cm**3*gauss/g')


#core_list = [21,70]
core_list = looper.get_all_nonzero()  #check what returns on get_all_nonzero
#frame_list=[125]
frame_list=[1,10,20,30,40,50,60,70,80,90,100,110,120,125]
fields=['density']
     
if 'this_looper' not in dir():
    directory = '/archive2/dcollins4096/Paper19/u05-r4-l4-128'
    this_looper = looper.core_looper(directory= directory,
                                     derived=[],
                                     sim_name = 'u05',
                                     out_prefix = 'pdf',
                                     target_frame = 125,
                                     frame_list = frame_list,
                                     core_list = core_list,
                                     fields_from_grid=['x','y','z']+fields
                                  )
    this_looper.get_target_indices(h5_name='datasets_small/u05_0125_peaklist.h5',
                                   bad_particle_list='datasets_small/bad_particles.h5')
    this_looper.get_tracks() 

import testing.early_mask as em
reload(em)

def toplot(prof,quan ='magnetic_field_strength'): 
    xbins = prof.x_bins
    bin_center = 0.5*(xbins[1:]+xbins[:-1]) 
    pdf = prof[quan]
    return xbins, bin_center,pdf

rm = rainbow_map(len(frame_list))

bbbA = []
bcenA = []
valsA = []

bbbC = []
bcenC = []
valsC = []

print("READY")
for index,frame in enumerate(frame_list):
    if 1: 
        ds = this_looper.load(frame=frame,derived=[em.add_tracer_density])
        em.add_tracer_density(ds)  
        ad = ds.all_data()
        # test new fields here with ad['field']
        deposit_tuple = ("deposit","target_particle_volume") 
        
    if 1:
        all_target_indices = np.concatenate([this_looper.target_indices[core_id] for core_id in core_list])
        ad.set_field_parameter('target_indices',all_target_indices)
        ad.set_field_parameter('mask_to_get',np.zeros_like(all_target_indices,dtype='int32'))
        #bins={'velocity_x':np.linspace(-25,25,64)}
        #bins['PotentialField']= np.linspace(-50,50,64)
        bins=None

        phase_all = yt.create_profile(ad,['density','magnetic_field_strength'],'cell_volume',weight_field=None, override_bins=bins)

        prof_alldata = yt.create_profile(ad,bin_fields=['density'],fields=['magnetic_field_strength'],weight_field='cell_volume', override_bins=bins)    
        prof_coresdata = yt.create_profile(ad,bin_fields=['density'],fields=['magnetic_field_strength'],weight_field='target_particle_volume', override_bins=bins) 
      
        bbb1, bcen1, vals1 = toplot(prof_alldata)
        bbb2, bcen2, vals2 = toplot(prof_coresdata)#,deposit_tuple[1])

        # ONE FRAME PDFS (one plot)
        if 0:
            bbbA.append(bbb1) 
            bcenA.append(bcen1)
            valsA.append(vals1)

            bbbC.append(bbb2)
            bcenC.append(bcen2)
            valsC.append(vals2)

        # 2 PLOTS PER FRAME (all_data, core_data)
        if 0: 
            outname = "plots_to_sort/Rho_B_%s.png"%frame
            fig,ax=plt.subplots(1,1)
            ax.plot(bcen1,vals1,c='k',linewidth=1.0)
            ax.plot(bcen2,vals2,c='g',linewidth=1.0)

            axbonk(ax,xlabel=r'$\rho$',ylabel=r'$\mid B \mid$',xscale='log',yscale='log')#,xlim=(1e-4,1e8),ylim=(1e-10,1e0))
            fig.savefig(outname)
            print(outname)
            plt.close(fig)

        # 4 PLOTS PER FRAME (all_data, core_data, phase-all_data + scatter plot) 
        if 1:
            p = yt.ParticlePhasePlot.from_profile(phase_all)
            p.set_xlim(1e-3,1e7)
            p.set_ylim(1e-1,1e4)  
            #p.set_cmap('cell_volume','arbre')
            p.set_zlim('cell_volume',1e-10,1e-1) 
            p.save()
            cell_plot = p.plots['cell_volume']
            this_axes = cell_plot.axes
            this_axes.plot(bcen1,vals1,c='k')
            this_axes.plot(bcen2,vals2,c='k',linestyle='dashed') 
 
            frame = index+2  #this only works if all frames are executed in order  
            brp.plot_particles(this_axes,frame) 
            p.save('4plots_%d'%frame)
           
# FOR ALL TIME PDFS (one plot)
if 0:
    fig,ax=plt.subplots(1,1)

    for i in range(len(frame_list)):
        ax.plot(bcenA[i],valsA[i],c=rm(i),linewidth=1.0)#,linestyle='dashed') 
        ax.plot(bcenC[i],valsC[i],c=rm(i),linewidth=0.4)#linestyle='dashed')
     
    outname = "plots_to_sort/PDFs.png"
    axbonk(ax,xlabel=r'$v$',ylabel=r'$V(v)$',xscale='log',yscale='log',xlim=(1e-4,1e8),ylim=(1e-10,1e0))
    #ax.set_title('%s'%frame)
    fig.savefig(outname)
    print(outname)
    plt.close(fig)

