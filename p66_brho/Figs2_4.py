
# <B> vs <n>
from starter2 import *
import data_locations as dl
from collections import defaultdict
import davetools
reload(davetools)

import testing.early_mask as em  #testing old method
reload(em)  #testing old method

import scipy
plt.close('all')


class magfield_density_tool(): 
    def __init__(self,this_looper):
        self.this_looper=this_looper
        
        self.mean_rho=defaultdict(list)
        self.mean_rho_py=defaultdict(list) 
        
        self.mean_field_comps=defaultdict(list) 
        self.mean_field_comps_py=defaultdict(list) 
        self.mean_field_x=defaultdict(list)
        self.mean_field_y=defaultdict(list)
        self.mean_field_z=defaultdict(list)

        self.mean_fieldOverRho = defaultdict(list) 
        self.mean_fieldOverRhoCheck = defaultdict(list) 
        self.mean_fieldOverRhoLogged = defaultdict(list)
          
        self.angle_mean = defaultdict(list)
        self.mean_cv=defaultdict(list)
        self.alphaSum = defaultdict(list)
        self.alphaProd = defaultdict(list)

        self.cores_used=[] 

        self.alpharr1 = np.empty([0],dtype=float)
        self.alpharr2 = np.empty([0],dtype=float)
        self.alpharr3 = np.empty([0],dtype=float)
        self.alpharr1_ad = np.empty([0],dtype=float)
        self.alpharr2_ad = np.empty([0],dtype=float)
        self.alpharr3_ad = np.empty([0],dtype=float)

        self.pearr1 = np.empty([0],dtype=float)
        self.pearr2 = np.empty([0],dtype=float)
        self.pearr3 = np.empty([0],dtype=float)

        
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
        self.core_list=core_list 
        
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            self.ms = ms 

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

                # <B(t=0)>
                bx_o = thtr.c([core_id],'magnetic_field_x')[mask,0]
                by_o = thtr.c([core_id],'magnetic_field_y')[mask,0]
                bz_o = thtr.c([core_id],'magnetic_field_z')[mask,0]
                bb_o = np.sqrt(bx_o*bx_o + by_o*by_o + bz_o*bz_o)
                
                BtdotBo = bx*bx_o + by*by_o + bz*bz_o 
                costheta = BtdotBo/(bb*bb_o)
                mean_cos_theta = (density*cell_volume*costheta).sum()/(density*cell_volume).sum()

                # EXPLORING THE ADDITIVE VS MULTIPLICATIVE INTEGRALS 
                bbsum = bb.sum()
                rhosum = density.sum()
                bblog = np.log(bb)
                rholog = np.log(density)

                fig,ax=plt.subplots(1,1)
                alpha1 = np.log(bbsum)/np.log(rhosum) 
                alpha2 = bblog.sum()/rholog.sum()
                        
                # EXPLORING OTHER B over RHO COMBINATIONS
                BRho = bb/density 
                density_log = np.log10(density)
                field_log = np.log10(bb)
                BRho_log = np.log10(BRho) 
                BovRho = density_log/field_log  #for lnB/lnRho

                # GETTING THE AVERAGES
                self.angle_mean[core_id].append(mean_cos_theta) 

                self.mean_field_comps[core_id].append((bb * cell_volume).sum()/cell_volume.sum())
                self.mean_field_x[core_id].append((bx * cell_volume).sum()/cell_volume.sum())
                self.mean_field_y[core_id].append((by * cell_volume).sum()/cell_volume.sum())
                self.mean_field_z[core_id].append((bz * cell_volume).sum()/cell_volume.sum())
                #self.mean_field_comps_py[core_id].append(bb.mean())  

                self.mean_rho[core_id].append((density * cell_volume).sum()/(cell_volume.sum()))  
                #self.mean_rho_py[core_id].append(density.mean())  
                self.mean_cv[core_id].append((cell_volume * cell_volume).sum()/(cell_volume.sum()))  #does this make sense to do...

                self.mean_fieldOverRho[core_id].append((BRho * cell_volume).sum()/(cell_volume.sum()))
                self.mean_fieldOverRhoCheck[core_id].append((BRho_log * cell_volume).sum()/(cell_volume.sum()))
                self.mean_fieldOverRhoLogged[core_id].append((BovRho * cell_volume).sum()/(cell_volume.sum()))

                self.alphaSum[core_id].append(alpha1)
                self.alphaProd[core_id].append(alpha2)


    def profiles(self,name=None):
            thtr = self.this_looper.tr 
            all_cores = np.unique(thtr.core_ids)
            core_list = all_cores

            for nf,theframe in enumerate(thtr.frames): 
                print('name ',name)
                fig,ax=plt.subplots(1,1)

                ds = self.this_looper.load(frame=theframe) 
                ad = ds.all_data()

                #all_target_indices = np.concatenate([self.this_looper.target_indices[core_id] for core_id in core_list])
                #ad.set_field_parameter('target_indices',all_target_indices)
                #ad.set_field_parameter('mask_to_get',np.zeros_like(all_target_indices,dtype='int32'))  

                if 0:  #ALL DATA
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
                    if name == 'u601':
                        self.alpharr1_ad = np.append(self.alpharr1_ad,alpha)
                    if name == 'u602':
                        self.alpharr2_ad = np.append(self.alpharr2_ad,alpha)
                    if name == 'u603':
                        self.alpharr3_ad = np.append(self.alpharr3_ad,alpha) 
                    Bavg_o = pfit[1]

                    X = np.linspace(the_x.min(),the_x.max(),num=len(the_x))  #short: -2, +3    
                    XX = 10 ** X
                    Y = 10 ** (pfit[0]*X + pfit[1])                 
                    #ax.plot(XX,Y,c='b',linewidth=1.0)
                    #ax.set_title(r'$\alpha = %.3f$'%alpha) 
                    

                if 1:  #ALL CORE DATA 
                    # the YT way:
                    #prof_cores = yt.create_profile(ad,bin_fields=['density'],fields=['magnetic_field_strength'],weight_field=deposit_tuple, override_bins=None)
                    #xbins = prof_cores.x_bins
                    #bin_center = 0.5*(xbins[1:]+xbins[:-1])
                    #pdf = prof_cores['magnetic_field_strength']
                    #the_xx = bin_center
                    #the_yy = pdf
                    #ax.plot(the_xx,the_yy,c='k',linewidth=1.0,linestyle='dashed')

                    # fields:
                    Cellvolume_all_cores_all_time = thtr.track_dict['cell_volume']  # shape: (particles,14) for u601
                    CellVolume = Cellvolume_all_cores_all_time[:,nf]  # shape: (particles,) for u601 

                    Rho_all_cores_all_time = thtr.track_dict['density'] 
                    Density = Rho_all_cores_all_time[:,nf]

                    Magfield_all_cores_all_time = thtr.track_dict['magnetic_field_strength']
                    Magfield = Magfield_all_cores_all_time[:,nf]
              
                    #BovRho_all_cores_all_time = Magfield_all_cores_all_time/Rho_all_cores_all_time     
                    #BovRho = BovRho_all_cores_all_time[:,nf]
                    BovRho = Magfield/Density  #SIMPLY

                    # 2D histograms:
                    cellvolume_bins = np.geomspace(1e-11,1e-6) 
                    density_bins = np.geomspace(1e-3,1e8)  #pdfs on overleaf
                    magfield_bins = np.geomspace(1e-1,1e4)
                     

# maybe there's a better way to do this, but perhaps this is the way to weight all our fields with cell volume
# these pdfs would correspond to alpha fit then average since it involves the particles of all cores at the same time
# (the black dot which essentially corresponds to the circled blue dot)
                    #theArray, binsX, binsY = np.histogram2d(Density, CellVolume, bins=[density_bins,cellvolume_bins],density=True)
                    theArray, binsX = np.histogram(BovRho,bins=density_bins, weights=CellVolume, density =True)
                     
                    # BinsX are the edges, center it by averaging:
                    bin_centers = 0.5*(binsX[1:]+binsX[:-1])
                    bin_widths = binsX[1:]-binsX[:-1]
                    #why do we need these next two lines again?
                    bin_centers.shape = bin_centers.size, 1   
                    bin_widths.shape = bin_widths.size, 1 

                    the_cxx = bin_centers  
                    the_cyy = theArray
                    
                    # attempts for 2D histograms:
                    #the_cyy = theArray.sum(axis=1) 
                    #the_cyy = (bin_centers * bin_widths * theArray).sum(axis=1) 
                    #the_cyy = (bin_widths * theArray).sum(axis=1) 

                    # it should make sense:
                    #whatsthis = the_cxx * the_cyy * bin_widths
                    #print('one or not',whatsthis)
                    #isitone = (the_cyy * bin_widths).sum()   #THINK!! or ask.. :/
                    #print('one or not',isitone)

                    ax.plot(the_cxx,the_cyy,c='k',linewidth=1.0)#, linestyle='dashed') 
                    #ax.set_xlim(1e-3,1e7)  #write in axbonk
                    ax.set_ylim(1e-10,1e1)  #write in axbonk 
                   
                    if 0:
                        # FOR HEAT MAP PURPOSES
                        #fig2,ax2 = plt.subplots(1,1) 
                        xbins = 0.5*(binsX[1:]+binsX[:-1])
                        ybins = 0.5*(binsY[1:]+binsY[:-1])
                        nx = len(xbins) ; ny=len(ybins)
                        # 2d array of xbins
                        TheX = np.r_[(ny)*[xbins]].transpose()
                        TheY = np.r_[(nx)*[ybins]]
                      

                        # all other heat maps have been for time in the x axis, something like this has been done:
                        # but here have our theArray for one time frame...and we want all time frames, let's see...

                        #hist = np.zeros([xbins.size,ybins.size])
                        #for ntime, time in enumerate(these_times):
                        #    thishist,bins = np.histogram(masses[ntime,:],bins=mass_bins_edge)
                        #    hist[ntime,:]=thishist
                        #print('len_hist ',hist.shape)  #49,49
                        #print('len_hist[nf,:] ',hist[nf,:].shape) #49
                        #print('len_theArray ',theArray.shape)  #49,49

                        cmap = copy.copy(mpl.cm.get_cmap("viridis"))
                        cmap.set_under('w')
                        minmin = theArray[theArray>0].min()
                        norm = mpl.colors.LogNorm(vmin=minmin,vmax=theArray.max())
                        ploot=ax.pcolormesh(TheX, TheY, theArray, cmap=cmap,norm=norm,shading='nearest')  # interpolation ~ shading

                    axbonk(ax,xscale='log',yscale='log') 
                    #ax.imshow(theArray,origin='lower',interpolation='nearest')  # imshow ~ plotting in Python
                    outname = 'PDF_BovRho_%s_%s'%(theframe,name)  
                    fig.savefig(outname)
                    print('save ',outname)
                    plt.close(fig)
                    
                #axbonk(ax,xlabel=r'$\rho/\rho_{o}$',ylabel=r'$|B|(\mu g)$',xscale='log',yscale='log') 
                #outname = 'Brho_profs_%s_%s'%(theframe,name) 
                #fig.savefig(outname)
                #print(outname)
                #plt.close(fig)



#import three_loopers_mountain_top as TLM
#import three_loopers_tenfour as TLTF
import three_loopers_six as TL6
if 'clobber' not in dir():
    clobber=True

if 'mag_den1' not in dir() or clobber:
    mag_den1=magfield_density_tool(TL6.loops['u601'])
    simname1 = 'u601'
    mag_den1.run()
    #mag_den1.profiles(simname1)

if 'mag_den2' not in dir() or clobber:
    mag_den2=magfield_density_tool(TL6.loops['u602'])
    simname2 = 'u602'
    mag_den2.run()
    #mag_den2.profiles(simname2)

if 'mag_den3' not in dir() or clobber:
    mag_den3=magfield_density_tool(TL6.loops['u603'])
    simname3 = 'u603'
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
        #print('tffs!',these_times)  #check the TIME  

        rhos = np.zeros([ntimes,ncores]) 
        fields = np.zeros([ntimes,ncores])
        fields_x = np.zeros([ntimes,ncores])
        fields_y = np.zeros([ntimes,ncores])
        fields_z = np.zeros([ntimes,ncores])
        fieldOvers = np.zeros([ntimes,ncores])
        fieldOversCheck = np.zeros([ntimes,ncores])
        fieldOversLog = np.zeros([ntimes,ncores])

        angles = np.zeros([ntimes,ncores])
        cvs = np.zeros([ntimes,ncores])
        alphaS = np.zeros([ntimes,ncores])
        alphaP = np.zeros([ntimes,ncores])


        # OPEN UP ALL THE FIGURES AND AXIS 
        # SINGLE PANEL 
        if 1:
            fig, ax1=plt.subplots(1,1) 
            fig0, ax0=plt.subplots(1,1)

        # OR MULTIPLE PANEL
        if 0:
            fig = plt.figure()
            #fig.text(0.365,0.03,r'$\rho/\rho_{o}$')

            ax1 = plt.subplot(331)
            ax2 = plt.subplot(332)
            ax3 = plt.subplot(333)  #ADDED 
            ax4 = plt.subplot(334)
            ax5 = plt.subplot(335)
            ax6 = plt.subplot(336)  #ADDED 
            ax7 = plt.subplot(337)
            ax8 = plt.subplot(338)
            ax9 = plt.subplot(339)  #ADDED 

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


        if 0:  # IF DOING <B> vs <n> BY FRAME, SINGLE PANEL
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
       

        the_v = np.empty([0],dtype=float)
        the_w = np.empty([0],dtype=float)

        the_x = np.empty([0],dtype=float)
        the_xX = np.empty([0],dtype=float)
        the_sxx = np.empty([0],dtype=float)
        the_a = np.empty([0],dtype=float)
        
        the_y = np.empty([0],dtype=float)
        the_syy = np.empty([0],dtype=float)
        the_b = np.empty([0],dtype=float)

        the_yx = np.empty([0],dtype=float)
        the_syxx = np.empty([0],dtype=float)
        the_c = np.empty([0],dtype=float)

        the_yl = np.empty([0],dtype=float)
        the_syll = np.empty([0],dtype=float)
        the_d = np.empty([0],dtype=float)

        the_yz = np.empty([0],dtype=float)
        the_yzs = np.empty([0],dtype=float)
        the_syzz = np.empty([0],dtype=float)
        the_syzzs = np.empty([0],dtype=float)
        the_e = np.empty([0],dtype=float)
      
        the_sxyz = np.empty([0],dtype=float) 
        the_xyz = np.empty([0],dtype=float) 

        the_z = np.empty([0],dtype=float) 
        the_cv = np.empty([0],dtype=float) 
 
        print('new sim')
        # MAKE THE FIELDS INTO A 2D ARRAY WE CAN PLOT
        for ncore,core_id in enumerate(tool.cores_used):
            this_rho = tool.mean_rho[core_id] 
            this_field = tool.mean_field_comps[core_id] 
            this_field_x = tool.mean_field_x[core_id]
            this_field_y = tool.mean_field_y[core_id]
            this_field_z = tool.mean_field_z[core_id]

            this_BRho = tool.mean_fieldOverRho[core_id]
            this_BRhoCheck = tool.mean_fieldOverRhoCheck[core_id]
            this_BRhoLog = tool.mean_fieldOverRhoLogged[core_id] 
            this_ang = tool.angle_mean[core_id]
            this_cv = tool.mean_cv[core_id]
             
            this_alphaS = tool.alphaSum[core_id]
            this_alphaP = tool.alphaProd[core_id]
            
            # passing the tool to the initiated field of zeros
            rhos[:,ncore]= this_rho 
            fields[:,ncore]= this_field
            fields_x[:,ncore]=this_field_x
            fields_y[:,ncore]=this_field_y
            fields_z[:,ncore]=this_field_z
            fieldOvers[:,ncore]= this_BRho
            fieldOversCheck[:,ncore]= this_BRhoCheck
            fieldOversLog[:,ncore] = this_BRhoLog
            angles[:,ncore]= np.arccos(this_ang)*180/np.pi
            cvs[:,ncore]= this_cv
            alphaS[:,ncore] = this_alphaS
            alphaP[:,ncore] = this_alphaP

            # passing it back to the original name
            this_alphaS = alphaS[:,ncore]
            this_alphaP = alphaP[:,ncore]

            this_rho = rhos[:,ncore] 
            the_xx = np.log10(this_rho)
            the_aa = this_rho
            the_x= np.append(the_x,the_xx)  
            #print('the_xX0',len(the_x))
            the_xX= np.append(the_xX,the_xx)
            #print('the_xX1',len(the_xX))
            the_xX= np.append(the_xX,the_xx)  
            #print('the_xX2',len(the_xX))
            the_xX= np.append(the_xX,the_xx)  #I know..but it works..  
            #print('the_xX3',len(the_xX))
            the_a= np.append(the_a,the_aa) 
           
            this_field = fields[:,ncore]
            the_yy = np.log10(this_field)
            the_bb = this_field
            the_y= np.append(the_y,the_yy)  
            the_b= np.append(the_b,the_bb)  
            
            this_field_x = abs(fields_x[:,ncore])  
            the_yxx = np.log10(this_field_x)
            the_cc = this_field_x
            the_yx =np.append(the_yx,the_yxx) 
            the_c = np.append(the_c,the_cc)
            
            this_field_y = abs(fields_y[:,ncore])  
            the_yll = np.log10(this_field_y)
            the_dd = this_field_y
            the_yl =np.append(the_yl,the_yll) 
            the_d = np.append(the_d,the_dd)

            this_field_z = abs(fields_z[:,ncore])  
            this_field_zz = fields_z[:,ncore]  
            the_yzz = np.log10(this_field_z)
            the_yzzs = np.log10(this_field_zz)
            the_ee = this_field_z
            the_yz =np.append(the_yz,the_yzz) 
            the_yzs =np.append(the_yzs,the_yzzs) 

            the_e = np.append(the_e,the_ee)
            the_xyz = np.append(the_xyz,the_yxx)
            the_xyz = np.append(the_xyz,the_yll) 
            the_xyz = np.append(the_xyz,the_yzz) 
          
            the_fieldOvRho = this_field/this_rho
            the_zz = np.log10(the_fieldOvRho)
            the_z = np.append(the_z,the_zz)
            
            this_BRho = fieldOvers[:,ncore]
            the_vv = np.log10(this_BRho)
            the_v= np.append(the_w,the_vv)  
            this_BRhoCheck = fieldOversCheck[:,ncore]
            this_BRhoLog = fieldOversLog[:,ncore]

            this_ang=angles[:,ncore]
            # TRY
            this_cv = cvs[:,ncore]
            the_cv=np.append(the_cv,this_cv)
      

            if 1: # FOR ALL TIME
                tmap = rainbow_map(len(this_rho))  #used to be[:-1]
                ctr = [tmap(n) for n in range(len(this_rho))]  #used to be[:-1]
                #ax1.scatter(this_field, this_field_z,c=ctr,marker='*')  
                #this_field = this_field_z
                #ax1.plot(this_field,this_field_z,linestyle='dashed')
                ax1.scatter(this_rho, this_field,c=ctr,marker='*',alpha=0.5)  
                #ax1.scatter(this_rho[0], this_field[0],c='b',marker='*')  
                #ax1.scatter(this_rho[-1], this_field[-1],c='r',marker='*')  
                #ax1.plot(this_rho,this_field,c=[0.5]*4)
                #ax1.set_xlabel(r'$B_avg$')
                #ax1.set_ylabel(r'$B_avg$')
                #ax1.set_xscale('log')
                #ax1.set_yscale('log')
                if 0: # IF NOT PLOTTING POWER LAW
                    outname_all='BzBTracks_%s'%simnames[nt]
                    fig.savefig(outname_all)
                    print("saved")

            if 0: # FOR ONE FRAME PER TIME; SINGLE PANEL 
                #xlims = 0.4,1.8
                #ylims = 0.0,5.0
                for i in range(len(axplts)): 
                    if 1:
                        if nt == 0:
                            color = 'r'
                        if nt == 1:
                            color = 'b'
                        if nt == 2:
                            color = 'g'
                    if 0:
                        axplts[i].scatter(this_alphaS[i],this_alphaP[i],c=color,marker='*')
                        outname_frame='AlphaSP_%s_%d'%(simnames[nt],i)
                        magfield_density_tool.labelled(axplts[i],xscale=None,yscale=None,xlabel='Sum',ylabel='Product',\
                                                       title=None, xlim=xlims,ylim=ylims)
                    if 1:
                        #axplts[i].scatter(this_field[i],this_field_x[i],c=color,marker='*',alpha=0.4)
                        #axplts[i].scatter(this_rho[i],this_field[i],c=color,marker='*',alpha=0.4)
                        axplts[i].scatter(this_rho[i],this_field_x[i],c='k',marker='x',alpha=0.4)
                        axplts[i].scatter(this_rho[i],this_field_y[i],c='k',marker='1',alpha=0.4)
                        axplts[i].scatter(this_rho[i],this_field_z[i],c='k',marker='*',alpha=0.4)
                        #axplts[i].scatter(this_rho[i],the_zz[i],c='g',alpha=0.2)

                        # IF PLOTTING THE POWER LAW, we can comment this out
                        #magfield_density_tool.labelled(axplts[i],xscale='log',yscale='log',xlabel=r'$\rho$',ylabel=r'$Bz$',\
                        #                               title=None, xlim=xlims,ylim=ylims)
                        #outname_frame='Scatter_LogBzvsRho_%s_%d'%(simnames[nt],i)
                    #figs[i].savefig(outname_frame)
                    #print("saved")
                

            if 0:  #for all extents I should combine these next ones with the extents in main definition
                   #maybe include under the next switch
                rho_extents=davetools.extents()
                rho_extents(the_a)
                magfield_extents = davetools.extents()
                magfield_extents(the_b)
                magfieldz_extents = davetools.extents()
                magfieldz_extents(the_c)

            if 0: # FOR ONE FRAME PER TIME; MULTIPLE PANEL ..in progress: I think it will be best if I end up moving this to the next if statement.
                tmap2 = rainbow_map(len(this_rho))
                c2 = [tmap2(n) for n in range(len(this_rho))]  
                for i in range(len(this_rho)):
                    if i in frames: 
                        lplots[i].scatter(this_rho[i],this_field[i],c=c2[i],marker='*')  #C2 gives me lots of notes, but it works... 
                        #lplots[n_time].plot(xx2,yy,c='k',linewidth=1.0)   #SEE BELOW
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

        print('max this_field_z',max(the_c))
        print('min this_field_z',min(the_c))
        if 0:
            plt.close(figs[i])

        the_w = the_y/the_x

        salpha_r = np.empty([0],dtype=float)
        salpha_xr = np.empty([0],dtype=float)
        salpha_yr = np.empty([0],dtype=float)
        salpha_zr = np.empty([0],dtype=float) 
        if 1:
            # PLOT THE POWER LAW: per frame 
            #numcores = len(the_x)/ncores 
            numcores = len(the_xX)/(ncores*3)   #for all comps 
            coreint = int(numcores)
            for i in range(coreint): 
                the_sx = the_x[i::coreint]
                the_sX = the_xX[i::coreint]  #for alpha of all comps...
                the_sy = the_y[i::coreint]
                the_syx = the_yx[i::coreint]
                the_syl = the_yl[i::coreint]
                the_syz = the_yz[i::coreint]
                the_sxyz = the_xyz[i::coreint]  #for alpha of all comps
                the_syzs = the_yzs[i::coreint]  #BZ SIGNED
                the_sz = the_z[i::coreint] 
                the_sw = the_w[i::coreint] 
                the_cvs = the_cv[i::coreint]

                # for power law for all time minus the last frame 
                #minusone = coreint-1
                #if i < minusone:  #for ALL frames 
                the_sxx= np.append(the_sxx,the_sx)   
                the_syy= np.append(the_syy,the_sy) 
                the_syzz= np.append(the_syzz,the_syz)
                the_syzzs = np.append(the_syzzs,the_syzs)  

                #sX = np.linspace(the_sx.min(),the_sx.max(),num=len(the_sx))  #short: -2, +3   
                sX = np.linspace(the_sX.min(),the_sX.max(),num=len(the_sX))  #short: -2, +3   
                #sX = np.linspace(the_sy.min(),the_sy.max(),num=len(the_sy))  #short: -2, +3   

                #spfit = np.polyfit(the_sy,the_syx,1)
                #spfit = np.polyfit(the_sy,the_syl,1)
                #spfit = np.polyfit(the_sy,the_syz,1)
                spfit = np.polyfit(the_sX,the_sxyz,1)  #all comps at once attempt
                #spfit = np.polyfit(the_sx,the_sy,1)
                spfit_x = np.polyfit(the_sx,the_syx,1)
                spfit_y = np.polyfit(the_sx,the_syl,1)
                spfit_z = np.polyfit(the_sx,the_syz,1)
                salpha = spfit[0]
                salpha_x = spfit_x[0]
                salpha_y = spfit_y[0]
                salpha_z = spfit_z[0]
                salpha_r = np.append(salpha_r,salpha)
                salpha_xr = np.append(salpha_xr,salpha_x)
                salpha_yr = np.append(salpha_yr,salpha_y)
                salpha_zr = np.append(salpha_zr,salpha_z)

                # pearsonR
                if 1:
                    xs = np.std(the_sx)
                    #ys = np.std(the_sy)
                    ys = np.std(the_syz)
                    if xs != 0 and ys != 0:
                        #pearX,pearY = scipy.stats.pearsonr(the_sx,the_sy)
                        pearX,pearY = scipy.stats.pearsonr(the_sx,the_syz)
                        #print('pearX',pearX)
                    else:
                        print("A zero encountered!!",xs,ys)
                
                if nt == 0:
                    mag_den1.alpharr1 = np.append(mag_den1.alpharr1,salpha) 
                    #print('601 bz alpha',mag_den1.alpharr1)
                    #mag_den1.pearr1 = np.append(mag_den1.pearr1,pearX)    
                if nt == 1:
                    mag_den2.alpharr2 = np.append(mag_den2.alpharr2,salpha) 
                    #print('602 bz alpha',mag_den2.alpharr2)
                    #mag_den2.pearr2 = np.append(mag_den2.pearr2,pearX)    
                if nt == 2:
                    mag_den3.alpharr3 = np.append(mag_den3.alpharr3,salpha) 
                    #print('603 bz alpha',mag_den3.alpharr3)
                    #mag_den3.pearr3 = np.append(mag_den3.pearr3,pearX)    
               
                sXX = 10 ** sX 
                sY = 10 ** (spfit[0]*sX + spfit[1])                
                sYx = 10 ** (spfit_x[0]*sX + spfit_x[1])                
                sYy = 10 ** (spfit_y[0]*sX + spfit_y[1])                
                sYz = 10 ** (spfit_z[0]*sX + spfit_z[1])                
                
                # PER PANEL
                # need x,y limits
                if 0:
                    xlims = 1e-1,1e8
                    ylims = 1e-2,1e4
                    print('inside plotting the power law')
                    print('i',i)
                    print('numcores',numcores)
                    print('len(the_xX)',len(the_xX))
                    axplts[i].plot(sXX,sY,c='k',linewidth=1.0)
                    #axplts[i].plot(sXX,sY,c='k',linewidth=1.0)
                    #axplts[i].plot(sXX,sYx,c='r',linestyle='dotted',linewidth=1.0)
                    #axplts[i].plot(sXX,sYy,c='r',linestyle='dashdot',linewidth=1.0)
                    #axplts[i].plot(sXX,sYz,c='r',linestyle='dashed',linewidth=1.0)

                    #pdb.set_trace()  #EDITTTT   
                    magfield_density_tool.labelled(axplts[i],xscale='log',yscale='log',xlabel=r'$<\rho>$',ylabel=r'$<Bx,By,Bz>$',\
                             xlim=xlims, ylim=ylims,title=r'$\alpha = %.3f$'%salpha)
                    #magfield_density_tool.labelled(axplts[i],xscale='log',yscale='log',xlabel=r'$<\rho>$',ylabel=r'$<B>$',\
                    #         xlim=xlims, ylim=ylims,title=r'$\alpha = %.3f, \alpha_x = %.3f, \alpha_y = %.3f, \alpha_z = %.3f$'%(salpha,salpha_x,salpha_y,salpha_z)) 
               
                    #outname_frame='Scatter_LogBxvsB_%s_%d'%(simnames[nt],i)
                    outname_frame='Scatter_LogBxyz3dvsRho_%s_%d'%(simnames[nt],i)
                    #outname_frame='Scatter_LogBzvsRho_%s_%d'%(simnames[nt],i)
                    figs[i].savefig(outname_frame)
                    print("saved ",i)
                    plt.close(figs[i])


                # MULTIPLE PANELS:
                if 0: 
                    if i in frames:
                        lplots[i].plot(sXX,sY,c='k',linewidth=1.0)

                # FOR HISTOGRAMS
                if 0:  
                    if nt == 0:
                        color = 'r'
                    if nt == 1:
                        color = 'b'
                    if nt == 2:
                        color = 'g'

                    # the_zmin -3.3661703693691116            the_zmax 1.0872747818424593
                    # the_wmin -73.79359150405679             the_wmax 82.79924734574196

                    the_bins = np.linspace(-10,10)
                    the_weights = None  #the_cvs 

                    # THE LN(Y/X)
                    the_lnArray, xbins = np.histogram(the_sz,bins=the_bins,weights=the_weights,density=True)
                    bin_lncenters = 0.5*(xbins[1:]+xbins[:-1])
                    the_lnX= bin_lncenters
                    the_lnY= the_lnArray
                    # THE LNY/LNX
                    the_lnlnArray, xxbins = np.histogram(the_sw,bins=the_bins,weights=the_weights,density=True)
                    bin_lnlncenters = 0.5*(xxbins[1:]+xxbins[:-1])
                    the_lnlnX= bin_lnlncenters
                    the_lnlnY= the_lnlnArray

                    axplts[i].plot(the_lnX,the_lnY,c=color,linewidth=1.0)
                    axplts[i].plot(the_lnlnX,the_lnlnY,c=color,linewidth=1.0,linestyle='dashed')
                    
                    #axplts[i].hist(the_sz, 50, density=False, histtype='step', color=color)  #EDIT color to reflect the small dcc py file
                    #axplts[i].hist(the_szz, 50, density=False, histtype='step', color=color,ls='--')  #EDIT color to reflect the small dcc py file
 
                    xlims = None
                    ylims = None  #check for limits
                    magfield_density_tool.labelled(axplts[i],xscale=None,yscale=None,xlabel=r'$log B/ \rho$, $logB/log\rho$',ylabel=r'$PDF$',\
                                       title=None, xlim=xlims, ylim=ylims)
                    # ensuring percentage for the y axis 
                    y_vals = axplts[i].get_yticks()
                    axplts[i].set_yticklabels(['{:.3f}'.format(x/len(the_lnX)) for x in y_vals])

                    outname_frame='TwoHistograms_AvgQty%s_%d'%(simnames[nt],i)
                    figs[i].savefig(outname_frame)
                    print("saved")
             

        if 0: 
            # SAVE THE MULTIPLE PANELS WITH THEIR LEAST SQUARE FIT
            outname = 'brhotffpanels_avgs_%s'%(simnames[nt])
            plt.savefig(outname)
            print("saved")

        if 0:
            alphaFile = open("alphaRecords.txt",'a')   
            #alphaFile.write("Sim %s Pears %f \n"%(simnames[nt],mag_den1.pearr1))  #how to write an array?
            print(mag_den1.pearr1)
            print(mag_den2.pearr2)
            print(mag_den3.pearr3)
            #alphaFile.write("Sim %s Alpha %f \n"%(simnames[nt],salpha))
            alphaFile.close()

        if 1:
            # PLOT THE POWER LAW: all time 
            pfit = np.polyfit(the_sxx,the_syy,1)  #B VS RHO
            #pfit = np.polyfit(the_syy,the_syzz,1)  #BZsigned VS B
            alpha = pfit[0]
            print('alpha',alpha)
            Bavg_o = pfit[1]
            # AND THE PEARSON R
            allxs = np.std(the_sxx)
            allys = np.std(the_syy)
            if xs != 0 and ys != 0:
                pearX,pearY = scipy.stats.pearsonr(the_sx,the_sy)
            else:
                print("A zero encountered!!",xs,ys)
            print('pearX_%s'%simnames[nt])
            print(pearX)

            X = np.linspace(the_sxx.min(),the_sxx.max(),num=len(the_sxx))  #short: +2, -3   
            #X = np.linspace(the_syy.min(),the_syy.max(),num=len(the_syy))  #short: +2, -3   
            XX = 10 ** X
            Y = 10 ** (pfit[0]*X + pfit[1])                
            
            if 1:
                ax1.plot(XX,Y,c='m',linewidth=1.0)
                if 1: # ORGANIZE YOUR PLOTTINGS!
                    xlabels = r'$\left\langle \rho/\rho_{o} \right\rangle$'
                    ylabels = r'$\left\langle\mid B \mid\right\rangle (\mu G)$'
                    xlims = 1e-1,1e8
                    ylims = 1e0,1e4
                    #ax1.set_title(r'$m = %.3f$'%alpha)
                    magfield_density_tool.labelled(ax1,xscale='log',yscale='log',xlabel=xlabels,ylabel=ylabels,\
                             xlim=xlims, ylim=ylims)#,title=r'$\alpha = %.3f$'%alpha)

                outname_all='BnTracks_pl_final%s'%simnames[nt]
                #outname_all='BzBTracks_%s'%simnames[nt]
                fig.savefig(outname_all)
                print("saved")


        # !!FOR LATER... #a tool for colors (tools/colors.py)!!! also look at make_core_cmap
        # PLOT ALL ALPHA COMPS VS TIME
        if 0:  
            outname = 'alphaPerFrame'  #swap order with below if for allCores
            outname = 'alphaPerFrame_allData'  #and comment next one out
            outname = 'alphaRecords'
            outname = 'pearRecords'
            outname = 'alpha_xyz3d_%s'%simnames[nt]
            outname_o = 'alphaxyz_vs_3d_scatt_%s'%simnames[nt]
            if nt == 0:
                therange = np.arange(0,1,0.075)
                #ax0.plot(therange, mag_den1.alpharr1_ad, c='g')                
                #ax0.plot(therange, mag_den1.alpharr1, c='g',linestyle='dashed')                
                #ax0.plot(therange, mag_den1.pearr1, c='g',linestyle='dashed')                
                ax1.plot(therange, salpha_r, c='k')                  
                ax1.plot(therange, salpha_xr, c='g',linestyle='dashed')
                ax1.plot(therange, salpha_yr, c='b',linestyle='dashed')
                ax1.plot(therange, salpha_zr, c='r',linestyle='dashed')

                ax0.plot(salpha_r,salpha_xr, c='g')                  
                ax0.plot(salpha_r,salpha_yr, c='b')                  
                ax0.plot(salpha_r,salpha_zr, c='r')                  
                ax0.scatter(salpha_r[2],salpha_xr[2], c='g',alpha=0.4)                  
                ax0.scatter(salpha_r[3],salpha_xr[3], c='g',alpha=0.4)                  
                ax0.scatter(salpha_r[8],salpha_xr[8], c='g',alpha=0.4)                  
                ax0.scatter(salpha_r[2],salpha_yr[2], c='b',alpha=0.4)                  
                ax0.scatter(salpha_r[3],salpha_yr[3], c='b',alpha=0.4)                  
                ax0.scatter(salpha_r[8],salpha_yr[8], c='b',alpha=0.4)                  
                ax0.scatter(salpha_r[2],salpha_zr[2], c='r',alpha=0.4)                  
                ax0.scatter(salpha_r[3],salpha_zr[3], c='r',alpha=0.4)                  
                ax0.scatter(salpha_r[8],salpha_zr[8], c='r',alpha=0.4)                  
                #ax0.plot(salpha_xr,salpha_yr, c='k')                  

            if nt == 1:
                therange = np.arange(0,1,0.075) 
                salpha_r = np.append(salpha_r,[np.nan]*1)
                salpha_xr = np.append(salpha_xr,[np.nan]*1)
                salpha_yr = np.append(salpha_yr,[np.nan]*1)
                salpha_zr = np.append(salpha_zr,[np.nan]*1)
                #mag_den2.alpharr2 = np.append(mag_den2.alpharr2,[np.nan]*1) 
                #mag_den2.alpharr2_ad = np.append(mag_den2.alpharr2_ad,[np.nan]*1) 
                #mag_den2.pearr2 = np.append(mag_den2.pearr2,[np.nan]*1) 
                #ax0.plot(therange, mag_den2.alpharr2_ad, c='b')
                #ax0.plot(therange, mag_den2.alpharr2, c='b',linestyle='dashed')
                #ax0.plot(therange, mag_den2.pearr2, c='b',linestyle='dashed')
                ax1.plot(therange, salpha_r, c='k')
                ax1.plot(therange, salpha_xr, c='g',linestyle='dashed')
                ax1.plot(therange, salpha_yr, c='b',linestyle='dashed')
                ax1.plot(therange, salpha_zr, c='r',linestyle='dashed')

                ax0.plot(salpha_r,salpha_xr, c='g')                  
                ax0.plot(salpha_r,salpha_yr, c='b')                  
                ax0.plot(salpha_r,salpha_zr, c='r')                  
                ax0.scatter(salpha_r[3],salpha_xr[3], c='g',alpha=0.4)                  
                ax0.scatter(salpha_r[8],salpha_xr[8], c='g',alpha=0.4)                  
                ax0.scatter(salpha_r[3],salpha_yr[3], c='b',alpha=0.4)                  
                ax0.scatter(salpha_r[8],salpha_yr[8], c='b',alpha=0.4)                  
                ax0.scatter(salpha_r[3],salpha_zr[3], c='r',alpha=0.4)                  
                ax0.scatter(salpha_r[8],salpha_zr[8], c='r',alpha=0.4)                  
                #ax0.plot(salpha_xr,salpha_yr, c='k')                  

            if nt == 2:
                therange = np.arange(0,1,0.075) 
                salpha_r = np.append(salpha_r,[np.nan]*2)
                salpha_xr = np.append(salpha_xr,[np.nan]*2)
                salpha_yr = np.append(salpha_yr,[np.nan]*2)
                salpha_zr = np.append(salpha_zr,[np.nan]*2)
                #mag_den3.alpharr3 = np.append(mag_den3.alpharr3,[np.nan]*2) 
                #mag_den3.alpharr3_ad = np.append(mag_den3.alpharr3_ad,[np.nan]*2) 
                #mag_den3.pearr3 = np.append(mag_den3.pearr3,[np.nan]*2) 
                #ax0.plot(therange, mag_den3.alpharr3_ad, c='m')
                #ax0.plot(therange, mag_den3.alpharr3, c='m',linestyle='dashed')
                #ax0.plot(therange, mag_den3.pearr3, c='m',linestyle='dashed')
                ax1.plot(therange, salpha_r, c='k')
                ax1.plot(therange, salpha_xr, c='g',linestyle='dashed')
                ax1.plot(therange, salpha_yr, c='b',linestyle='dashed')
                ax1.plot(therange, salpha_zr, c='r',linestyle='dashed')

                ax0.plot(salpha_r,salpha_xr, c='g')                  
                ax0.plot(salpha_r,salpha_yr, c='b')                  
                ax0.plot(salpha_r,salpha_zr, c='r')                  
                ax0.scatter(salpha_r[2],salpha_xr[2], c='g',alpha=0.4)                  
                ax0.scatter(salpha_r[5],salpha_xr[5], c='g',alpha=0.4)                  
                ax0.scatter(salpha_r[8],salpha_xr[8], c='g',alpha=0.4)                  
                ax0.scatter(salpha_r[2],salpha_yr[2], c='b',alpha=0.4)                  
                ax0.scatter(salpha_r[5],salpha_yr[5], c='b',alpha=0.4)                  
                ax0.scatter(salpha_r[8],salpha_yr[8], c='b',alpha=0.4)                  
                ax0.scatter(salpha_r[2],salpha_zr[2], c='r',alpha=0.4)                  
                ax0.scatter(salpha_r[5],salpha_zr[5], c='r',alpha=0.4)                  
                ax0.scatter(salpha_r[8],salpha_zr[8], c='r',alpha=0.4)                  
                #ax0.plot(salpha_xr,salpha_yr, c='k')                  
                
            ax1.set_xlabel(r'$t/t_{\rm{ff}}$')
            #ax0.set_ylabel(r'$\alpha_{B\rho}$')
            #ax0.set_ylabel(r'$R_{B\rho}$')
            ax1.set_ylabel(r'$\alpha,\alpha_x,\alpha_y,\alpha_z$')
            #ax0.set_ylim(-.1,1.0)
            #mag_den3.labelled(ax0,xlabel=r'$t_{\rm{ff}}$',ylabel=r'$R_{B\rho}$',ylim=ylim,title=None) 
            fig.savefig(outname)
            print(outname," saved")
            plt.close(fig)

            #ax0.set_xlabel(r'$\alpha_x$')
            ax0.set_xlabel(r'$\alpha$')
            ax0.set_ylabel(r'$\alpha_x,\alpha_y,\alpha_z$')
            fig0.savefig(outname_o)
            print(outname_o," saved")
            plt.close(fig0)

    if 0:
        print("calling alphaFieldDensity.py")
        afd.axisforbox(ax0)  #plots the boxplot and returns

        fig0.savefig(outname)   #saves all three plots 
        print(outname," saved")
        plt.close(fig0)


