
# <B> vs <n>
from starter2 import *
import data_locations as dl
from collections import defaultdict
import davetools
reload(davetools)

#import testing.early_mask as em  #testing old method
#reload(em)  #testing old method

import scipy
plt.close('all')


class magfield_density_tool(): 
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[] 
        
        self.mean_rho=defaultdict(list)
        self.mean_field_comps=defaultdict(list) 
        
    def labelled(ax,xscale=None,yscale=None,xlabel=None,ylabel=None,\
                 xlim=None,ylim=None,title=None,linthreshx=0.1,linthreshy=0.1):  
        if xscale and yscale != None:
            ax.set_xscale(xscale)
            ax.set_yscale(yscale)
        if xscale == 'symlog':
            ax.set_xscale(xscale,linthreshx=linthreshx)
        if yscale == 'symlog':
            ax.set_yscale(yscale,linthreshy=linthreshy)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel) 
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_title(title)

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
        
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            self.ms = ms 

            self.cores_used.append(core_id)
            self.times = thtr.times
            asort =  np.argsort(self.times) 

            for nf,frame in enumerate(thtr.frames): 
                mask = ms.compute_unique_mask(core_id, dx=1./2048,frame=nf)

                density = thtr.c([core_id],'density')[mask,nf]
                cell_volume = thtr.c([core_id],'cell_volume')[mask,nf]
                cell_mass = density * cell_volume 

                bx = thtr.c([core_id],'magnetic_field_x')[mask,nf]
                by = thtr.c([core_id],'magnetic_field_y')[mask,nf]
                bz = thtr.c([core_id],'magnetic_field_z')[mask,nf]                
                bb = np.sqrt(bx*bx+by*by+bz*bz) 
 
                # GETTING THE AVERAGES
                self.mean_field_comps[core_id].append((bb * density * cell_volume).sum()/cell_mass.sum())  
                self.mean_rho[core_id].append((density * cell_volume).sum()/(cell_volume.sum()))  


import three_loopers_six as TL6
if 'clobber' not in dir():
    clobber=True

if 'mag_den1' not in dir() or clobber:
    mag_den1=magfield_density_tool(TL6.loops['u601'])
    simname1 = 'u601'
    mag_den1.run()

if 'mag_den2' not in dir() or clobber:
    mag_den2=magfield_density_tool(TL6.loops['u602'])
    simname2 = 'u602'
    mag_den2.run()

if 'mag_den3' not in dir() or clobber:
    mag_den3=magfield_density_tool(TL6.loops['u603'])
    simname3 = 'u603'
    mag_den3.run()
simnames = [simname1, simname2, simname3]

'''
TURN ON WHICH FIGURE?
'''
# IF FIGURE 2 
if 0:
    time = 'all_time'
# IF FIGURE 4
if 1:
    time = 'per_frame'
    meanbear = []


for nt,tool in enumerate([mag_den1,mag_den2,mag_den3]):   
    # SET UP THE VARIABLE
    G=1620./(4*np.pi)
    ncores = len(tool.cores_used)
    ntimes = tool.times.size

    rhos = np.zeros([ntimes,ncores]) 
    fields = np.zeros([ntimes,ncores])

    # OPEN UP ALL THE FIGURES AND AXIS 
    # SINGLE PANEL; Fig 2. 
    if time == 'all_time':
        fig, ax1=plt.subplots(1,1) 
        fig0, ax0=plt.subplots(1,1)

    # MULTIPLE PANEL; Fig 4.
    if time == 'per_frame':
        fig = plt.figure()

        ax1 = plt.subplot(331)
        ax2 = plt.subplot(332)
        ax3 = plt.subplot(333) 
        ax4 = plt.subplot(334)
        ax5 = plt.subplot(335)
        ax6 = plt.subplot(336) 
        ax7 = plt.subplot(337)
        ax8 = plt.subplot(338)
        ax9 = plt.subplot(339) 

        fig.subplots_adjust(wspace=0, hspace=0)

        if nt == 0:
            frames = [1,2,3,5,6,7,8,10,12]
            lplots = [0,ax1,ax2,ax3,0,ax4,ax5,ax6,ax7,0,ax8,0,ax9,0] 
        if nt == 1:
            frames = [1,2,3,4,5,6,7,9,11]
            lplots = [0,ax1,ax2,ax3,ax4,ax5,ax6,ax7,0,ax8,0,ax9,0] 
        if nt == 2:
            frames = [1,2,3,4,5,6,7,9,10]
            lplots = [0,ax1,ax2,ax3,ax4,ax5,ax6,ax7,0,ax8,ax9,0] 

    the_x = np.empty([0],dtype=float)
    the_sxx = np.empty([0],dtype=float)
    the_a = np.empty([0],dtype=float)
    
    the_y = np.empty([0],dtype=float)
    the_syy = np.empty([0],dtype=float)
    the_b = np.empty([0],dtype=float)

    # MAKE THE FIELDS INTO A 2D ARRAY WE CAN PLOT
    for ncore,core_id in enumerate(tool.cores_used):
        this_rho = tool.mean_rho[core_id] 
        this_field = tool.mean_field_comps[core_id] 

        # passing the tool to the initiated field of zeros
        rhos[:,ncore]= this_rho 
        fields[:,ncore]= this_field

        # passing it back to the original name...
        this_rho = rhos[:,ncore] 
        the_xx = np.log10(this_rho)
        the_x= np.append(the_x,the_xx)  
        the_aa = this_rho
        the_a= np.append(the_a,the_aa)
 
        this_field = fields[:,ncore]
        the_yy = np.log10(this_field)
        the_y= np.append(the_y,the_yy)  
        the_bb = this_field
        the_b= np.append(the_b,the_bb)
        
        if time == 'all_time': 
            tmap = rainbow_map(len(this_rho)) 
            ctr = [tmap(n) for n in range(len(this_rho))]  
            ax1.scatter(this_rho, this_field,c=ctr,alpha=0.5) #marker='*'  
            pdb.set_trace()
 
        if time == 'per_frame': 
            rho_extents=davetools.extents()
            rho_extents(the_a)
            magfield_extents = davetools.extents()
            magfield_extents(the_b)

            tmap2 = rainbow_map(len(this_rho))
            c2 = [tmap2(n) for n in range(len(this_rho))]  
            for i in range(len(this_rho)):
                if i in frames: 
                    lplots[i].scatter(this_rho[i],this_field[i],c=[c2[i]]*this_rho[i].size,s=2,marker='o')  #C2 gives me lots of notes, but it works... 
                    #pdb.set_trace()
                    magfield_density_tool.labelled(lplots[i],xscale='log',yscale='log',xlabel=None,ylabel=None,\
                                                   title=None, xlim=rho_extents.minmax,ylim=magfield_extents.minmax)
                    ax2.tick_params(axis='y',labelleft=False)
                    ax3.tick_params(axis='y',labelleft=False)
                    ax4.set_ylabel(r'$\left\langle\mid B \mid\right\rangle (\mu G)$')
                    ax5.tick_params(axis='y',labelleft=False)
                    ax6.tick_params(axis='y',labelleft=False)
                    ax8.tick_params(axis='y',labelleft=False)
                    ax8.set_xlabel(r'$\left\langle \rho/\rho_{o} \right\rangle$')
                    ax9.tick_params(axis='y',labelleft=False)


    numcores = len(the_x)/ncores 
    coreint = int(numcores)
    for i in range(coreint): 
        the_sx = the_x[i::coreint]
        the_sy = the_y[i::coreint]

        the_sxx= np.append(the_sxx,the_sx)   
        the_syy= np.append(the_syy,the_sy) 

        sX = np.linspace(the_sx.min(),the_sx.max(),num=len(the_sx))  #short: -2, +3   

        spfit = np.polyfit(the_sx,the_sy,1)
        salpha = spfit[0]
        meanbear.append(salpha)

        sXX = 10 ** sX 
        sY = 10 ** (spfit[0]*sX + spfit[1])                
       

        # PLOT THE POWER LAW: per frame 
        if time == 'per_frame': 
            if i in frames:
                lplots[i].plot(sXX,sY,c='grey',linewidth=0.8)

            outname = 'brhotffpanels_avgs_%s'%(simnames[nt])
            plt.savefig(outname)
            #plt.close(figs[i])
            print("saved")

            # THE PEARSON R 
            xs = np.std(the_sx)
            ys = np.std(the_sy)
            if xs != 0 and ys != 0:
                pearX,pearY = scipy.stats.pearsonr(the_sx,the_sy)
            else:
                print("A zero encountered!!",xs,ys)
    print('meanbear',meanbear)


    # PLOT THE POWER LAW: all time 
    if time == 'all_time':
        pfit = np.polyfit(the_sxx,the_syy,1) 
        alpha = pfit[0]
        Bavg_o = pfit[1]

        X = np.linspace(the_sxx.min(),the_sxx.max(),num=len(the_sxx))  #short: +2, -3   
        XX = 10 ** X
        Y = 10 ** (pfit[0]*X + pfit[1])                
        
        ax1.plot(XX,Y,c='grey',linewidth=2.0, linestyle='dashed')
        xlabels = r'$\left\langle \rho/\rho_{o} \right\rangle$'
        ylabels = r'$\left\langle\mid B \mid\right\rangle (\mu G)$'
        xlims = 1e-1,1e8
        ylims = 1e0,1e4
        magfield_density_tool.labelled(ax1,xscale='log',yscale='log',xlabel=xlabels,ylabel=ylabels,\
                 xlim=xlims, ylim=ylims)

        outname_all='BnTracks_pl_final%s'%simnames[nt]
        fig.savefig(outname_all)
        print("saved")
        
        # THE PEARSON R 
        allxs = np.std(the_sxx)
        allys = np.std(the_syy)
        if allxs != 0 and allys != 0:
            pearX,pearY = scipy.stats.pearsonr(the_sxx,the_syy)
        else:
            print("A zero encountered!!",xs,ys)


