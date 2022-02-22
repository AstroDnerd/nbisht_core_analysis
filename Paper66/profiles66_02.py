
'''
PART OF P19_NEWSCRIPTS REPOSITORY:
This script must be placed on: ~/p19_newscripts
With: starter1.py, starter2.py on same directory.
Started as: ~/p19_newscripts/examples/mask_profile.py
In conjuction with: Brho_particles66.py


A SCRIPT FOR PROFILES, PHASES, PLOT
'''

from starter2 import *
from starter1 import *
#import Brho_particles66 as brp 
import data_locations as dl


def _BoverRho(field,data):
    return data['magnetic_field_strength']/data['density']
yt.add_field('_BoverRho',function=_BoverRho,sampling_type='cell',units='cm**3*gauss/g')

def _OmegaOverRho(field,data):
    return data['vorticity_magnitude']/data['density']
yt.add_field('_OmegaOverRho',function=_OmegaOverRho,sampling_type='cell',units='cm**3/(g*s)')

def _BoverRhoVort(field,data):
    return data['magnetic_field_strength']/(data['density'] * data['vorticity_magnitude'])
yt.add_field('_BoverRhoVort',function=_BoverRhoVort,sampling_type='cell',units='cm**3*gauss*s/g')




if 'this_simname' not in dir():
    this_simname = 'u05'
    
if 'this_looper' not in dir():
    file_list=glob.glob(dl.every_ten[this_simname])
    this_looper=looper.core_looper(directory=dl.sims[this_simname])
    for nfile,fname in enumerate(file_list):
        this_looper.load_loop(fname)
        print( "Reading file %d of %d"%(nfile,len(file_list)))
    thtr = this_looper.tr
    thtr.sort_time()
    all_cores = np.unique(thtr.core_ids)

core_list = all_cores
#core_list = looper.get_all_nonzero(dl.n_particles[this_simname])  #CAREFUL, check looper.py !!!
frame_list = dl.frame_list[this_simname]
fields=['density']

# FOR TESTING
#core_list = [3,5] 
#frame_list = [1]


# WILL REPLACING THIS WITH THE ABOVE WORK?? YES, data_puller.py does the process below once and for all to save us time 
'''
if 'this_looper' not in dir():
    directory = dl.sims[this_simname]
    this_looper = looper.core_looper(directory= directory,
                                     derived=[],
                                     sim_name = this_simname,
                                     out_prefix = 'sort_plots/%s_test'%this_simname,
                                     target_frame = dl.target_frame_list[this_simname],
                                     frame_list = frame_list,
                                     core_list = core_list,
                                     fields_from_grid=['x','y','z']+fields
                                  )
    this_looper.get_target_indices(h5_name=dl.peak_list[this_simname],
                                   bad_particle_list=dl.bad_particles[this_simname])
    this_looper.get_tracks() 
'''


import testing.early_mask as em
reload(em)

def toplot(prof,quan='vorticity_magnitude'):
    xbins = prof.x_bins
    bin_center = 0.5*(xbins[1:]+xbins[:-1]) 
    pdf = prof[quan]
    return xbins, bin_center,pdf
# MAYBE ADD TO toplot
    #bin_widths = xbins[1:]-xbins[:-1]
    #pdf = pdf/bin_widths
    #return xbins, bin_center,pdf,bin_widths

rm = rainbow_map(len(frame_list))

bbbA = []
bcenA = []
valsA = []

bbbC = []
bcenC = []
valsC = []

print("READY")
for index,frame in enumerate(frame_list):
    if 1:  #skip 
        ds = this_looper.load(frame=frame,derived=[em.add_tracer_density])
        em.add_tracer_density(ds)  
        ad = ds.all_data()
        # test new fields here with ad['field']
        deposit_tuple = ("deposit","target_particle_volume") 
        
    if 1:
        all_target_indices = np.concatenate([this_looper.target_indices[core_id] for core_id in core_list])
        ad.set_field_parameter('target_indices',all_target_indices)
        ad.set_field_parameter('mask_to_get',np.zeros_like(all_target_indices,dtype='int32')) 
        bins=None

        #phase_all = yt.create_profile(ad,['density','magnetic_field_strength'],'cell_volume',weight_field=None, override_bins=bins)
        phase_all = yt.create_profile(ad,['density','vorticity_magnitude'],'cell_volume',weight_field=None, override_bins=bins)

        prof_alldata = yt.create_profile(ad,bin_fields=['density'],fields=['vorticity_magnitude'],weight_field='cell_volume', override_bins=bins)
        prof_coresdata = yt.create_profile(ad,bin_fields=['density'],fields=['vorticity_magnitude'],weight_field=deposit_tuple, override_bins=bins) 
      
        #prof_alldata = yt.create_profile(ad,bin_fields=['vorticity_magnitude'],fields=['cell_volume'],weight_field=None, override_bins=bins)
        #prof_coresdata = yt.create_profile(ad,bin_fields=['vorticity_magnitude'],fields=[deposit_tuple],weight_field=None, override_bins=bins) 

        bbb1, bcen1, vals1 = toplot(prof_alldata)
        bbb2, bcen2, vals2 = toplot(prof_coresdata)#, quan=deposit_tuple[1])

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
            outname = "BoverRho_Vort_%s_%s.png"%(frame,this_simname) 
            fig,ax=plt.subplots(1,1)
            ax.plot(bcen1,vals1,c='k',linewidth=1.0)
            ax.plot(bcen2,vals2,c='b',linewidth=1.0, linestyle='dashed')

            #axbonk(ax,xlabel=r'$B/\rho$',ylabel=r'$Vorticity$',xscale='log',yscale='log',xlim=(1e-3,1e4),ylim=(1e2,1e6))
            #axbonk(ax,xlabel=r'$\rho$',ylabel=r'$\mid B \mid$',xscale='log',yscale='log')#,xlim=(1e-4,1e8),ylim=(1e-10,1e0))
            axbonk(ax,xlabel=r'$\nabla x v$',ylabel=r'$V(\nabla x v)$',xscale='log',yscale='log',xlim=(1e-4,1e8),ylim=(1e-10,1e0)) 
            fig.savefig(outname)
            print(outname)
            plt.close(fig)

        # 4 PLOTS PER FRAME (all_data, core_data, phase-all_data + scatter plot) 
        if 1:
            p = yt.ParticlePhasePlot.from_profile(phase_all)   
            #p.set_ylim(1e0,1e5)  # vorticity VS B/rho 
            #p.set_ylim(1e-2,1e6)  # vorticity/rho VS B/rho 
            p.set_ylim(1e0,1e5)  # vorticity VS B,rho 
            
            #p.set_xlim(1e-4,1e4)  # vorticity,/rho VS B/rho 
            #p.set_xlim(1e-1,1e4)  # vorticity VS B
            p.set_xlim(1e-2,1e7)  # vorticity VS B,rho
            
            ##p.set_cmap('cell_volume','arbre') 
            p.set_zlim('cell_volume',1e-9,1e-2) # vorticity VS B/rho, rho
            #p.set_zlim('cell_volume',1e-10,1e-2)  # vorticity,/rho VS B,/rho   
            p.save()
            
            cell_plot = p.plots['cell_volume']
            this_axes = cell_plot.axes
            this_axes.plot(bcen1,vals1,c='k')
            this_axes.plot(bcen2,vals2,c='k',linestyle='dashed') 
             
            #nframe = index  #works only if all frames are in order; check what is best +1,+2,0...find the general trend...? 
            #brp.plot_particles(this_axes,nframe) 
            p.save('threeplots_%d_%s'%(frame,this_simname)) 

# FOR ALL TIME PDFS (one plot)
if 0:
    fig,ax=plt.subplots(1,1)

    for i in range(len(frame_list)):
        ax.plot(bcenA[i],valsA[i],c=rm(i),linewidth=1.0)#,linestyle='dashed') 
        ax.plot(bcenC[i],valsC[i],c=rm(i),linewidth=0.4)#linestyle='dashed')
     
    outname = "sort_plots/RhoPDFs_%s.png"%this_simname
    axbonk(ax,xlabel=r'$\rho/\rho_{o}$',ylabel=r'$V(\rho)$',xscale='log',yscale='log',xlim=(1e-4,1e8),ylim=(1e-10,1e0))
    #axbonk(ax,xlabel=r'$B$',ylabel=r'$V(B)$',xscale='log',yscale='log',xlim=(1e-4,1e8),ylim=(1e-10,1e0))
    #ax.set_title('%s'%frame)
    fig.savefig(outname)
    print(outname)
    plt.close(fig)

