
#
from starter2 import *
import data_locations as dl
from collections import defaultdict

import davetools
reload(davetools)

plt.close('all')


class magfield_density_tool():  # was mass_tool()
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.mean_rho=defaultdict(list)
        self.mean_rho_py=defaultdict(list)

        self.mean_field=defaultdict(list)
        self.mean_field_py=defaultdict(list)
        
        self.mean_field_comps=defaultdict(list) 
        self.mean_field_comps_py=defaultdict(list) 
          
        self.angle_mean = defaultdict(list)
        self.cores_used=[] 


    def run(self,core_list=None):
        dx=1./2048
        nx = 1./dx
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        print('core_list',len(core_list))
        thtr.sort_time()
        tsorted = thtr.times
        self.core_list=core_list  #initiated?
        
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            self.ms = ms  #initiated?

            #print('go ', core_id)
            self.cores_used.append(core_id)
            self.times = thtr.times


            for nf,frame in enumerate(thtr.frames): 
                mask = ms.compute_unique_mask(core_id, dx=1./2048,frame=nf)

                density = thtr.c([core_id],'density')[mask,nf]
                cell_volume = thtr.c([core_id],'cell_volume')[mask,nf]

                #magfield = thtr.c([core_id],'magnetic_field_stregnth')[mask,nf]  #ERROR, CHECK[ip]
                bx = thtr.c([core_id],'magnetic_field_x')[mask,nf]
                by = thtr.c([core_id],'magnetic_field_y')[mask,nf]
                bz = thtr.c([core_id],'magnetic_field_z')[mask,nf]                
                bb = np.sqrt(bx*bx+by*by+bz*bz) 

                # <B(t=0)> : CHECK[  ]
                bx_o = thtr.c([core_id],'magnetic_field_x')[mask,0]
                by_o = thtr.c([core_id],'magnetic_field_y')[mask,0]
                bz_o = thtr.c([core_id],'magnetic_field_z')[mask,0]
                bb_o = np.sqrt(bx_o*bx_o + by_o*by_o + bz_o*bz_o)
                
                BtdotBo = bx*bx_o + by*by_o + bz*bz_o 
                costheta = BtdotBo/(bb*bb_o)

                mean_cos_theta = (density*cell_volume*costheta).sum()/(density*cell_volume).sum()
               
                #self.mean_field[core_id].append((magfield * cell_volume).sum()/cell_volume.sum())
                #self.mean_field_py[core_id].append(magfield.mean())

                self.mean_field_comps[core_id].append((bb * cell_volume).sum()/cell_volume.sum())
                #self.mean_field_comps_py[core_id].append(bb_comps.mean())  

                self.mean_rho[core_id].append((density * cell_volume).sum()/(cell_volume.sum()))  
                #self.mean_rho_py[core_id].append(density.mean())  

                self.angle_mean[core_id].append(mean_cos_theta)
                #print('good') 


import three_loopers_mountain_top as TLM
if 'clobber' not in dir():
    clobber=True
if 'mag_den1' not in dir() or clobber:
    mag_den1=magfield_density_tool(TLM.loops['u301'])
    simname1 = 'u301'
    mag_den1.run()
if 'mag_den2' not in dir() or clobber:
    mag_den2=magfield_density_tool(TLM.loops['u302'])
    simname2 = 'u302'
    mag_den2.run()
if 'mag_den3' not in dir() or clobber:
    mag_den3=magfield_density_tool(TLM.loops['u303'])
    simname3 = 'u303'
    mag_den3.run()
simnames = [simname1, simname2, simname3]


#
# Heat Map with Sample Tracks
#
if 0:
    fig,ax=plt.subplots(1,3, figsize=(12,4))
    axes=ax.flatten()

if 1:
    for nt,tool in enumerate([mag_den1,mag_den2,mag_den3]):
        # SET UP THE VARIABLE
        G=1620./(4*np.pi)
        tff_global = np.sqrt(3*np.pi/(32*G*1))
        ncores = len(tool.cores_used)
        ntimes = tool.times.size
        these_times = tool.times/tff_global

        rhos = np.zeros([ntimes,ncores])
        fields = np.zeros([ntimes,ncores])
        angles = np.zeros([ntimes,ncores])
        outname='BoBang_tracks_%s'%simnames[nt]

        if 1: # just the tracks 
            fig, ax1=plt.subplots(1,1)

        # MAKE THE FIELDS INTO A 2D ARRAY WE CAN PLOT
        for ncore,core_id in enumerate(tool.cores_used):
            this_rho = tool.mean_rho[core_id] 
            this_field = tool.mean_field_comps[core_id] 
            this_ang = tool.angle_mean[core_id]
 
            rhos[:,ncore]= np.log10(this_rho)
            fields[:,ncore]= np.log10(this_field)
            angles[:,ncore]= np.arccos(this_ang)*180/np.pi

            this_rho = rhos[:,ncore] 
            this_field = fields[:,ncore]

            if 1: # (just) all the tracks
                this_ang=angles[:,ncore]
                ax1.plot(these_times,this_ang,c=[0.5]*4)  
        if 1:
            ax1.plot([0,1],[90,90],'--',c='k') 
            fig.savefig(outname)
            print('saved')
''' 
        # PLOT A FEW OF THE TRACKS
        nc = len(tool.cores_used)
        take_a_few = ((nc-1)*np.random.random(10)).astype('int')
        for ncore,core_id in enumerate(nar(tool.cores_used)[take_a_few]): 
            this_ang=angles[:,ncore]
            axes[nt].plot(these_times,this_ang,c=[0.5]*4)  


        # MAKE A 2D HISTOGRAM
        #There are probably better ways to do this.
        angle_bins_edge = np.linspace(0,180,32)
        angle_bins_edge = np.linspace(0,180,32)
        xbins = these_times
        ybins = 0.5*(angle_bins_edge[1:]+angle_bins_edge[:-1])
        nx = len(xbins) ; ny=len(ybins)
        TheX = np.r_[(ny)*[xbins]].transpose()
        TheY = np.r_[(nx)*[ybins]]

        hist = np.zeros([xbins.size,ybins.size])
        for ntime, time in enumerate(these_times):
            thishist,bins = np.histogram(angles[ntime,:],bins=angle_bins_edge)
            hist[ntime,:]=thishist

        # PLOT THE MAP
        # First set up the color map, and force it to be white in areas where theres no points
        cmap = copy.copy(mpl.cm.get_cmap("viridis"))
        cmap.set_under('w')
        minmin = hist[hist>0].min()

        # SET UP THE COLORBAR
        norm = mpl.colors.LogNorm(vmin=1,vmax=33)
        ploot=axes[nt].pcolormesh(TheX, TheY, hist, cmap=cmap,norm=norm,shading='nearest')
        axbonk(axes[nt],ylabel=None,xlabel=r'$t/t_{\rm{ff}}$',yscale='linear', ylim=[0,180])
    axes[0].set_ylabel(r'$<\theta>$')
    fig.colorbar(ploot)
    fig.savefig(outname)
    print(outname)
'''
