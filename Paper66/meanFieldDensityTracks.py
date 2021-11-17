
# <B> vs <n>
from starter2 import *
import data_locations as dl
from collections import defaultdict
import davetools
reload(davetools)

import testing.early_mask as em  #testing old method
reload(em)  #testing old method
plt.close('all')


class magfield_density_tool(): 
    def __init__(self,this_looper):
        self.this_looper=this_looper
        
        self.mean_rho=defaultdict(list)
        self.mean_rho_py=defaultdict(list) 
        
        self.mean_field_comps=defaultdict(list) 
        self.mean_field_comps_py=defaultdict(list) 
          
        self.angle_mean = defaultdict(list)
        self.cores_used=[] 

        self.alpharr1 = np.empty([0],dtype=float)
        self.alpharr2 = np.empty([0],dtype=float)
        self.alpharr3 = np.empty([0],dtype=float)
        self.alpharr1_ad = np.empty([0],dtype=float)
        self.alpharr2_ad = np.empty([0],dtype=float)
        self.alpharr3_ad = np.empty([0],dtype=float)

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
        self.core_list=core_list  #initiated?
        
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            self.ms = ms  #initiated?

            #print('go ', core_id)
            self.cores_used.append(core_id)
            self.times = thtr.times
            asort =  np.argsort(self.times)  #added

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
                self.angle_mean[core_id].append(mean_cos_theta)
                
                self.mean_field_comps[core_id].append((bb * cell_volume).sum()/cell_volume.sum())
                #self.mean_field_comps_py[core_id].append(bb.mean())  

                self.mean_rho[core_id].append((density * cell_volume).sum()/(cell_volume.sum()))  
                #self.mean_rho_py[core_id].append(density.mean())  

    def profiles(self,name=None):
            thtr = self.this_looper.tr 
            all_cores = np.unique(thtr.core_ids)
            core_list = all_cores

            for nf,theframe in enumerate(thtr.frames): 
              
                fig,ax=plt.subplots(1,1)

                ds = self.this_looper.load(frame=theframe)  #derived=[em.add_tracer_density]), was testing old method
                #em.add_tracer_density(ds)
                ad = ds.all_data()
                #deposit_tuple = ("deposit","target_particle_volume") 

                #all_target_indices = np.concatenate([self.this_looper.target_indices[core_id] for core_id in core_list])
                #ad.set_field_parameter('target_indices',all_target_indices)
                #ad.set_field_parameter('mask_to_get',np.zeros_like(all_target_indices,dtype='int32'))  

                if 1:  #ALL DATA
                    prof_ad = yt.create_profile(ad,bin_fields=['density'],fields=['magnetic_field_strength'],weight_field='cell_volume', override_bins=None)

                    xbins = prof_ad.x_bins
                    bin_center = 0.5*(xbins[1:]+xbins[:-1])
                    pdf = prof_ad['magnetic_field_strength']

                    the_xx = bin_center
                    the_x = np.log10(ad['density'].v) 
                    the_yy = pdf
                    the_y = np.log10(ad['magnetic_field_strength'].v)

                    ax.plot(the_xx,the_yy,c='k',linewidth=1.0)
                    
                    # PLOT THE POWER LAW 
                    pfit = np.polyfit(the_x,the_y,1)
                    alpha = pfit[0]
                    if name == 'u401':
                        self.alpharr1_ad = np.append(self.alpharr1_ad,alpha)
                    if name == 'u402':
                        self.alpharr2_ad = np.append(self.alpharr2_ad,alpha)
                    if name == 'u403':
                        self.alpharr3_ad = np.append(self.alpharr3_ad,alpha) 
                    Bavg_o = pfit[1]

                    X = np.linspace(the_x.min(),the_x.max(),num=len(the_x))  #short: -2, +3    
                    XX = 10 ** X
                    Y = 10 ** (pfit[0]*X + pfit[1])                 
                    #ax.plot(XX,Y,c='b',linewidth=1.0)
                    #ax.set_title(r'$\alpha = %.3f$'%alpha) 
                    

                if 1:  #ALL CORE DATA
                 
                    #prof_cores = yt.create_profile(ad,bin_fields=['density'],fields=['magnetic_field_strength'],weight_field=deposit_tuple, override_bins=None)
                    #xbins = prof_cores.x_bins
                    #bin_center = 0.5*(xbins[1:]+xbins[:-1])
                    #pdf = prof_cores['magnetic_field_strength']
                    #the_xx = bin_center
                    #the_yy = pdf
                    #ax.plot(the_xx,the_yy,c='k',linewidth=1.0,linestyle='dashed')

                    CellVolume_all_cores_all_time = thtr.track_dict['cell_volume'] 
                    CellVolume = CellVolume_all_cores_all_time[:,nf]

                    Rho_all_cores_all_time = thtr.track_dict['density'] 
                    Density = Rho_all_cores_all_time[:,nf]

                    Magfield_all_cores_all_time = thtr.track_dict['magnetic_field_strength']
                    Magfield = Magfield_all_cores_all_time[:,nf]
                    
                    density_bins = np.geomspace(1e-3,1e7)
                    magfield_bins = np.geomspace(1e-1,1e4)
                    theArray, binsX, binsY = np.histogram2d(Density, Magfield, bins=[density_bins,magfield_bins], density=True) 
               
                    # BinsX are the edges, center it by averaging
                    bin_centers = 0.5*(binsX[1:]+binsX[:-1])
                    bin_widths = binsX[1:]-binsX[:-1]
                    bin_centers.shape = bin_centers.size, 1 
                    bin_widths.shape = bin_widths.size, 1 

                    the_cxx = bin_centers  
                    #the_cyy = theArray.mean(axis=1)  #TRY .SUM(axis1)
                    BRho = (bin_centers * bin_widths * theArray).sum(axis=1)          
                    ax.plot(the_cxx,BRho,c='k',linewidth=1.0, linestyle='dashed') 
                    ax.set_xlim(1e-3,1e7)  #write in axbonk
                    ax.set_ylim(1e-2,1e4)  #write in axbonk 
                   

                    # REVIEW THE FOLLOWING CODE UP TO # *
                    #fig2,ax2 = plt.subplots(1,1) 
                    #xbins = 0.5*(binsX[1:]+binsX[:-1])
                    #ybins = 0.5*(binsY[1:]+binsY[:-1])
                    #nx = len(xbins) ; ny=len(ybins)
                    # 2d array of xbins
                    #TheX = np.r_[(ny)*[xbins]].transpose()
                    #TheY = np.r_[(nx)*[ybins]]

                    #cmap = copy.copy(mpl.cm.get_cmap("viridis"))
                    #cmap.set_under('w')
                    #minmin = theArray[theArray>0].min()
                    #norm = mpl.colors.LogNorm(vmin=minmin,vmax=theArray.max())
                    #ploot=ax2.pcolormesh(TheX, TheY, theArray, cmap=cmap,norm=norm,shading='nearest')  # interpolation ~ shading
                    #ax2.set_yscale('log')
                    #ax2.set_xscale('log')
                    #ax2.imshow(theArray,origin='lower',interpolation='nearest')  # imshow ~ plotting in Python
                    #fig2.savefig('BRho_order')
                    # *
                    
                axbonk(ax,xlabel=r'$\rho/\rho_{o}$',ylabel=r'$|B|(\mu g)$',xscale='log',yscale='log') 
                outname = 'Brho_profs_%s_%s'%(theframe,name) 
                fig.savefig(outname)
                print(outname)
                plt.close(fig)
          


#import three_loopers_mountain_top as TLM
import three_loopers_tenfour as TLTF
if 'clobber' not in dir():
    clobber=True

if 'mag_den1' not in dir() or clobber:
    mag_den1=magfield_density_tool(TLTF.loops['u401'])
    simname1 = 'u401'
    mag_den1.run()
    #mag_den1.profiles(simname1)

if 'mag_den2' not in dir() or clobber:
    mag_den2=magfield_density_tool(TLTF.loops['u402'])
    simname2 = 'u402'
    mag_den2.run()
    #mag_den2.profiles(simname2)

if 'mag_den3' not in dir() or clobber:
    mag_den3=magfield_density_tool(TLTF.loops['u403'])
    simname3 = 'u403'
    mag_den3.run()
    #mag_den3.profiles(simname3)
simnames = [simname1, simname2, simname3]


# UNCOMMENT IF OVERLAYING ALL PLOTS
#import alphaFieldDensity as afd
if 1:
    # TO PLOT FIGUERS OF ALL THREE SIMS AT ONCE
    #fig0, ax0=plt.subplots(1,1) 

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


        # OPEN UP ALL THE FIGURES AND AXIS
        fig, ax1=plt.subplots(1,1) 
        if 0:  # if doing <B> vs <n> by frame
            fig2,ax2 = plt.subplots(1,1) 
            fig3,ax3 = plt.subplots(1,1) 
            fig4,ax4 = plt.subplots(1,1) 
            fig5,ax5 = plt.subplots(1,1) 
            fig6,ax6 = plt.subplots(1,1) 
            fig7,ax7 = plt.subplots(1,1) 
            fig8,ax8 = plt.subplots(1,1) 
            fig9,ax9 = plt.subplots(1,1) 
            fig10,ax10 = plt.subplots(1,1) 
            fig11,ax11 = plt.subplots(1,1) 
            fig12,ax12 = plt.subplots(1,1) 
            fig13,ax13 = plt.subplots(1,1)
            fig14,ax14 = plt.subplots(1,1)
            fig15,ax15 = plt.subplots(1,1)
            
            if nt == 0:
                axplts = [ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12,ax13,ax14,ax15]
                figs = [fig2,fig3,fig4,fig5,fig6,fig7,fig8,fig9,fig10,fig11,fig12,fig13,fig14,fig15]
            if nt == 1: 
                axplts = [ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12,ax13,ax14]
                figs = [fig2,fig3,fig4,fig5,fig6,fig7,fig8,fig9,fig10,fig11,fig12,fig13,fig14]
            if nt == 2: 
                axplts = [ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12,ax13]
                figs = [fig2,fig3,fig4,fig5,fig6,fig7,fig8,fig9,fig10,fig11,fig12,fig13]
        
        the_x = np.empty([0],dtype=float)
        the_y = np.empty([0],dtype=float)
        
        # MAKE THE FIELDS INTO A 2D ARRAY WE CAN PLOT
        for ncore,core_id in enumerate(tool.cores_used):
            this_rho = tool.mean_rho[core_id] 
            this_field = tool.mean_field_comps[core_id] 
            this_ang = tool.angle_mean[core_id]
 
            rhos[:,ncore]= this_rho  # was np.log10()
            fields[:,ncore]= this_field
            angles[:,ncore]= np.arccos(this_ang)*180/np.pi

            this_rho = rhos[:,ncore] 
            the_xx = np.log10(this_rho)
            the_x= np.append(the_x,the_xx) 
           
            this_field = fields[:,ncore]
            the_yy = np.log10(this_field)
            the_y= np.append(the_y,the_yy)  
           
            this_ang=angles[:,ncore]
            
            if 1: # for all frames
                tmap = rainbow_map(len(this_rho))
                ctr = [tmap(n) for n in range(len(this_rho))]
                ax1.scatter(this_rho, this_field,c=ctr,marker='*')  
                #ax1.scatter(this_rho[0], this_field[0],c='b',marker='*')  
                #ax1.scatter(this_rho[-1], this_field[-1],c='r',marker='*')  
                #ax1.plot(this_rho,this_field,c=[0.5]*4)

            if 0: 
                for i in range(len(axplts)):
                    axplts[i].scatter(this_rho[i],this_field[i],c='g',marker='*')

        xlims = 10e-2,10e7
        ylims = 10e-1,10e3
        if 0:
            # PLOT THE POWER LAW: per frame 
            numcores = len(the_x)/ncores 
            coreint = int(numcores)
            for i in range(coreint): 
                the_sx = the_x[i::coreint]
                the_sy = the_y[i::coreint]
                sX = np.linspace(the_sx.min(),the_sx.max(),num=len(the_sx))  #short: -2, +3   

                spfit = np.polyfit(the_sx,the_sy,1)
                salpha = spfit[0]
                sBavg_o = spfit[1]
                if nt == 0:
                    mag_den1.alpharr1 = np.append(mag_den1.alpharr1,salpha) 
                if nt == 1:
                    mag_den2.alpharr2 = np.append(mag_den2.alpharr2,salpha) 
                if nt == 2:
                    mag_den3.alpharr3 = np.append(mag_den3.alpharr3,salpha) 
               
                sXX = 10 ** sX 
                sY = 10 ** (spfit[0]*sX + spfit[1])                
 
                #axplts[i].plot(sXX,sY,c='k',linewidth=1.0)
                #magfield_density_tool.labelled(axplts[i],xscale='log',yscale='log',xlabel=r'$<\rho>$',ylabel=r'$<B>$',\
                #         xlim=xlims, ylim=ylims,title=r'$\alpha = %.3f$'%salpha)

                outname_frame='BnFrameTracks_pl_%s_%d'%(simnames[nt],i)
                #figs[i].savefig(outname_frame)
                print("saved ",i)
                plt.close('all')
 
                if 0:
                    alphaFile = open("alphaRecords.txt",'a')   
                    alphaFile.write("Sim %s Alpha %f \n"%(simnames[nt],salpha))
                    alphaFile.close()

        if 1:
            # PLOT THE POWER LAW: all frames 
            pfit = np.polyfit(the_x,the_y,1)
            alpha = pfit[0]
            Bavg_o = pfit[1]

            X = np.linspace(the_x.min()+2,the_x.max()-3,num=len(the_x))  #short: -2, +3   
            XX = 10 ** X
            Y = 10 ** (pfit[0]*X + pfit[1])                
                                                           
            ax1.plot(XX,Y,c='k',linewidth=1.0)
            magfield_density_tool.labelled(ax1,xscale='log',yscale='log',xlabel=r'$<\rho>$',ylabel=r'$<B>$',\
                     xlim=xlims, ylim=ylims)#,title=r'$\alpha = %.3f$'%alpha)

            outname_all='BnTracks_plxs_%s'%simnames[nt]
            fig.savefig(outname_all)
            print("saved")


        # !!FOR LATER... #a tool for colors (tools/colors.py)!!! also look at make_core_cmap
        if 0:  
            outname = 'alphaPerFrame'  #swap order with below if for allCores
            outname = 'alphaPerFrame_allData'  #and comment next one out
            outname = 'alphaRecords'
            if nt == 0:
                therange = np.arange(0,1,0.075)
                ax0.plot(therange, mag_den1.alpharr1_ad, c='g')                
                ax0.plot(therange, mag_den1.alpharr1, c='g',linestyle='dashed')                
            if nt == 1:
                therange = np.arange(0,1,0.075)
                #therange = np.arange(0.075,1,0.075)
                mag_den2.alpharr2 = np.append(mag_den2.alpharr2,[np.nan]*1) 
                mag_den2.alpharr2_ad = np.append(mag_den2.alpharr2_ad,[np.nan]*1) 
                ax0.plot(therange, mag_den2.alpharr2_ad, c='b')
                ax0.plot(therange, mag_den2.alpharr2, c='b',linestyle='dashed')
            if nt == 2:
                therange = np.arange(0,1,0.075)
                #therange = np.arange(0.15,1,0.075)
                mag_den3.alpharr3 = np.append(mag_den3.alpharr3,[np.nan]*2) 
                mag_den3.alpharr3_ad = np.append(mag_den3.alpharr3_ad,[np.nan]*2) 
                ax0.plot(therange, mag_den3.alpharr3_ad, c='m')
                ax0.plot(therange, mag_den3.alpharr3, c='m',linestyle='dashed')
                
                ax0.set_xlabel(r'$t_{\rm{ff}}$')
                ax0.set_ylabel(r'$\alpha_{B\rho}$')
                #mag_den3.labelled(ax0,xlabel=r'$t_{\rm{ff}}$',ylabel=r'$\alpha_B$',ylim=None,title=None) 

    if 0:
        print("calling alphaFieldDensity.py")
        afd.axisforbox(ax0)  #plots the boxplot and returns

        fig0.savefig(outname)   #saves all three plots 
        print(outname," saved")
        plt.close(fig0)


