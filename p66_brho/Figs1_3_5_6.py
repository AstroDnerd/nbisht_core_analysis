
# FOR THREE PLOTS OF PAPER66
'''
 PART OF P19_NEWSCRIPTS REPOSITORY:
This script must be placed on: ~/p19_newscripts
With: starter1.py, starter2.py on same directory.
 notes: 
 for debug purposes, long_list = long_list[:3]

    B = B0 rho^beta --> lnB = ln(B_0*rho^beta) -->
    lnB = beta*ln(rho) + lnB_0
'''

# - - - - - - - - - - - - - - - - - - - - 
from starter2 import *
import data_locations as dl
import davetools
reload(davetools)

import scipy
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick

import math
plt.close('all')


# - - - - - - - - - - - - - - - - - - - - 
# THE CLASS
class BRho_tool():
    def __init__(self,this_looper): 
        self.this_looper = this_looper
        self.betarr = np.empty([0],dtype=float)
        self.pearsonr = np.empty([0],dtype=float)

        self.pear0 = np.empty([0],dtype=float)
        self.pear1 = np.empty([0],dtype=float)
        self.pear2 = np.empty([0],dtype=float)
        self.pear3 = np.empty([0],dtype=float)
        self.pear4 = np.empty([0],dtype=float)
        self.pear5 = np.empty([0],dtype=float)
        self.pear6 = np.empty([0],dtype=float)
        self.pear7 = np.empty([0],dtype=float)
        self.pear8 = np.empty([0],dtype=float)
        self.pear9 = np.empty([0],dtype=float)
        self.pear10 = np.empty([0],dtype=float)
        self.pear11 = np.empty([0],dtype=float)
        self.pear12 = np.empty([0],dtype=float)
        self.pear13 = np.empty([0],dtype=float)
        self.pear14 = np.empty([0],dtype=float)

        self.bear0 = np.empty([0],dtype=float)
        self.bear1 = np.empty([0],dtype=float)
        self.bear2 = np.empty([0],dtype=float)
        self.bear3 = np.empty([0],dtype=float)
        self.bear4 = np.empty([0],dtype=float)
        self.bear5 = np.empty([0],dtype=float)
        self.bear6 = np.empty([0],dtype=float)
        self.bear7 = np.empty([0],dtype=float)
        self.bear8 = np.empty([0],dtype=float)
        self.bear9 = np.empty([0],dtype=float)
        self.bear10 = np.empty([0],dtype=float)
        self.bear11 = np.empty([0],dtype=float)
        self.bear12 = np.empty([0],dtype=float)
        self.bear13 = np.empty([0],dtype=float)
        #self.bear14 = np.empty([0],dtype=float)

        self.dear0 = np.empty([0],dtype=float)
        self.dear1 = np.empty([0],dtype=float)
        self.dear2 = np.empty([0],dtype=float)
        self.dear3 = np.empty([0],dtype=float)
        self.dear4 = np.empty([0],dtype=float)
        self.dear5 = np.empty([0],dtype=float)
        self.dear6 = np.empty([0],dtype=float)
        self.dear7 = np.empty([0],dtype=float)
        self.dear8 = np.empty([0],dtype=float)
        self.dear9 = np.empty([0],dtype=float)
        self.dear10 = np.empty([0],dtype=float)
        self.dear11 = np.empty([0],dtype=float)
        self.dear12 = np.empty([0],dtype=float)
        self.dear13 = np.empty([0],dtype=float)

        self.brms = np.empty([0],dtype=float)
        self.vrms = np.empty([0],dtype=float)
        self.vrmsq_bmean = np.empty([0],dtype=float)

        self.pears = [] 
        self.bears = []
        self.dears = []
   
        self.bears = [self.bear0,self.bear1,self.bear2,self.bear3,self.bear4,self.bear5,self.bear6,\
                      self.bear7,self.bear8,self.bear9,self.bear10,self.bear11,self.bear12,self.bear13]
        self.pears = [self.pear0,self.pear1,self.pear2,self.pear3,self.pear4,self.pear5,self.pear6,\
                      self.pear7,self.pear8,self.pear9,self.pear10,self.pear11,self.pear12,self.pear13]
        self.dears = [self.dear0,self.dear1,self.dear2,self.dear3,self.dear4,self.dear5,self.dear6,\
                      self.dear7,self.dear8,self.dear9,self.dear10,self.dear11,self.dear12,self.dear13]

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


    def run(self,name,thtr,timing,typeplot,tfflabs,
            fig=None,ax1=None,ax2=None,ax3=None,ax4=None,ax5=None,ax7=None,ax8=None,
            lplots=None,core_list=None,core_val=None,simframes=None):

        # FIG. 6
        def pearsonR(the_x, the_y):
            # LOGGED
            # NOTE: the "else" statements are a holder for the pearsonR's that are nans;
            # these need to be accounted for rightfully. suggestion: reverse the time,core loops
            xs = np.std(the_x)  
            ys = np.std(the_y)

            # SCIPY PEARSONr:
            if xs != 0 and ys != 0:
                pearX,pearY = scipy.stats.pearsonr(the_x,the_y)   
            else:
                print("A zero encountered!!",xs,ys)
                pearX = 0 
            self.pearsonr = np.append(self.pearsonr,pearX)  
       
        # FIG 1, 3, 5, 6. 
        def scatterplots(fig=None,ax1=None,ax2=None,ax3=None,ax4=None,ax5=None,ax7=None,ax8=None,
                         lplots=None,xx2=None,yy=None,yym=None):
         
            # FIG 1.
            if typeplot == 'scatter_plot' and core_id == core_val:
                if n_time <= asort[-1]:
                    ax1.scatter(density[mask,n_time],magfield[mask,n_time],c=c,label=thtr.times[n_time],s=0.1)

                xlims = 10e-3,10e6      
                ylims = 10e-2,10e3      
                if n_time == asort[-1]: 
                    ax1.plot(xx2,yy,c='k',linewidth=1.0)
                    ax1.plot(xx2,yym,'--',c='grey',linewidth=1.0)
                    BRho_tool.labelled(ax1,xscale='log',yscale='log',xlabel=r'$\rho/\rho_{o}$',ylabel=r'$\mid B \mid (\mu G)$',
                                       xlim=xlims, ylim=ylims)
            # FIG 3. 
            if typeplot == 'frame_scatters':
                if n_time in simframes:   
                    lplots[n_time].scatter(density[mask,n_time],magfield[mask,n_time],c=c,label=thtr.times[n_time],s=0.1)  #edit thtr.times          
                    lplots[n_time].plot(xx2,yy,c='k',linewidth=1.0) #c=c[0] for colors
                    print("scatters plotted")             
                  
                    #print('nTIME',n_time)
                    BRho_tool.labelled(lplots[n_time],xscale='log',yscale='log',xlabel=None,ylabel=None, 
                                       xlim=rho_extents.minmax, ylim=magfield_extents.minmax) 
                    ax2.tick_params(axis='y',labelleft=False)
                    ax4.set_ylabel(r'$\mid B \mid (\mu G)$')
                    ax5.tick_params(axis='y',labelleft=False)
                    ax8.tick_params(axis='y',labelleft=False) 

                if n_time == asort[-1]:                
                    tmap2 = rainbow_map(len(self.betarr)) 
                    c2 = [tmap2(n) for n in range(len(self.betarr[:-1]))] 

                    if name == 'u601': 
                        the_range = np.arange(0,1,0.075) 
                        ax3.scatter(the_range[:-1],self.betarr[:-1],c=c2)   #[2:],
                        ax3.plot(the_range[:-1],self.betarr[:-1],c='k',linewidth=1.0)  
                        #ax3.set_title(r"$\beta = 0.2, \alpha_{cf}$")
                        ax3.plot([0,1],[0.20,0.20],'--',c=[0.5]*4)
                    if name == 'u602':
                        the_range = np.arange(0.075,1,0.075) 
                        ax3.scatter(the_range[:-1],self.betarr[:-1],c=c2)  
                        ax3.plot(the_range[:-1],self.betarr[:-1],c='k',linewidth=1.0)  
                        #ax3.set_title(r"$\beta = 2.0, \alpha_{cf}$")
                        ax3.plot([0,1],[0.39,0.39],'--',c=[0.5]*4)
                    if name == 'u603':
                        the_range = np.arange(0.075,0.9,0.075) 
                        ax3.scatter(the_range[:-1],self.betarr[:-1],c=c2)  
                        ax3.plot(the_range[:-1],self.betarr[:-1],c='k',linewidth=1.0)   
                        #ax3.set_title(r"$\beta = 20, \alpha_{cf}$")
                        ax3.plot([0,1],[0.42,0.42],'--',c=[0.5]*4)
                    ax3.plot([0,1],[0,0],c=[0.5]*4) 
                 
                    tff_mod = [0.0,0.0,0.3,0.6,0.9]  #this should match multiplelocator
                    ax3.set_xticklabels(tff_mod)   
                    ax3.xaxis.set_major_locator(plt.MultipleLocator(0.3))

                    ax3.yaxis.tick_right()
                    ax3.set_xlabel(r'$t/t_{\rm{ff}}$')  #\rm{  },to avoid italization 
                    ax3.set_yscale('linear')
                    ax3.set_ylim(-0.5,0.5)
          
                    self.betarr = np.empty([0],dtype=float)
   
                    outname = 'brhotffpanels_c%04d_%s'%(core_id,name)  
                    plt.savefig(outname)  # CAREFUL WITH FIG VS PLT
                    print("saved "+outname) 

            # FIG. 5
            if typeplot == 'box_plot' and n_time == asort[-1]:   
                self.bear0 = np.append(self.bear0,self.betarr[0])
                self.bear1 = np.append(self.bear1,self.betarr[1])
                self.bear2 = np.append(self.bear2,self.betarr[2]) 
                self.bear3 = np.append(self.bear3,self.betarr[3]) 
                self.bear4 = np.append(self.bear4,self.betarr[4]) 
                self.bear5 = np.append(self.bear5,self.betarr[5]) 
                self.bear6 = np.append(self.bear6,self.betarr[6])  
                self.bear7 = np.append(self.bear7,self.betarr[7]) 
                self.bear8 = np.append(self.bear8,self.betarr[8]) 
                self.bear9 = np.append(self.bear9,self.betarr[9])
                self.bear10 = np.append(self.bear10,self.betarr[10]) 
                self.bear11 = np.append(self.bear11,self.betarr[11])
                if name =='u601' or name =='u602': 
                    self.bear12 = np.append(self.bear12,self.betarr[12])    
                if name == 'u601':
                    self.bear13 = np.append(self.bear13,self.betarr[13]) 

                if 1:
                    if core_id == core_list[-1]:
                        if name == 'u601':
                            meanbear = [0.207432,0.201813,0.234106,0.230126,0.160165,0.074837,0.107986,0.079870,0.030422,0.043188,0.101638,0.176109,0.353048,0.456282]    
                        if name == 'u602':
                            meanbear = [0.273494,0.265940,0.329641,0.385815,0.344329,0.284348,0.203154,0.173078,0.160007,0.166019,0.240440,0.277389,0.479462]
                        if name == 'u603': 
                            meanbear = [0.279555,0.336840,0.355424,0.325673,0.268460,0.232444,0.258084,0.252883,0.206200,0.239367,0.283474,0.617871]
                        self.bear0 = np.append(self.bear0,(meanbear[0]))
                        self.bear1 = np.append(self.bear1,(meanbear[1]))
                        self.bear2 = np.append(self.bear2,(meanbear[2]))
                        self.bear3 = np.append(self.bear3,(meanbear[3]))
                        self.bear4 = np.append(self.bear4,(meanbear[4]))
                        self.bear5 = np.append(self.bear5,(meanbear[5]))
                        self.bear6 = np.append(self.bear6,(meanbear[6]))
                        self.bear7 = np.append(self.bear7,(meanbear[7]))
                        self.bear8 = np.append(self.bear8,(meanbear[8]))
                        self.bear9 = np.append(self.bear9,(meanbear[9]))
                        self.bear10 = np.append(self.bear10,(meanbear[10]))
                        self.bear11 = np.append(self.bear11,(meanbear[11]))
                        if name =='u601' or name =='u602': 
                            self.bear12 = np.append(self.bear12,(meanbear[12]))
                        if name == 'u601':
                            self.bear13 = np.append(self.bear13,(meanbear[13]))

                self.bears = [self.bear0,self.bear1,self.bear2,self.bear3,self.bear4,self.bear5,self.bear6,\
                              self.bear7,self.bear8,self.bear9,self.bear10,self.bear11,self.bear12,self.bear13]
                self.betarr = np.empty([0],dtype=float)  


            # FIG 6.
            if typeplot == 'vio_plot' and n_time == asort[-1]:     
                self.pear0 = np.append(self.pear0,self.pearsonr[0])
                self.pear1 = np.append(self.pear1,self.pearsonr[1]) 
                self.pear2 = np.append(self.pear2,self.pearsonr[2]) 
                self.pear3 = np.append(self.pear3,self.pearsonr[3]) 
                self.pear4 = np.append(self.pear4,self.pearsonr[4]) 
                self.pear5 = np.append(self.pear5,self.pearsonr[5]) 
                self.pear6 = np.append(self.pear6,self.pearsonr[6])  
                self.pear7 = np.append(self.pear7,self.pearsonr[7]) 
                self.pear8 = np.append(self.pear8,self.pearsonr[8]) 
                self.pear9 = np.append(self.pear9,self.pearsonr[9]) 
                self.pear10 = np.append(self.pear10,self.pearsonr[10]) 
                self.pear11 = np.append(self.pear11,self.pearsonr[11]) 
                if name == 'u601' or name == 'u602':
                    self.pear12 = np.append(self.pear12,self.pearsonr[12])
                if name == 'u601':
                    self.pear13 = np.append(self.pear13,self.pearsonr[13])  

                if 1:
                    if core_id == core_list[-1]:
                        if name == 'u601':
                            meanpear = [0.46115629,0.386104,0.5580705,0.53245971,0.40396877,0.2024436,0.28630761,0.24313097,0.09931512,\
                                        0.17075329, 0.38030015, 0.57047563,0.78074059, 0.77001102]
                        if name == 'u602':
                            meanpear = [0.63225779, 0.59682206, 0.72420179, 0.80303703, 0.74460867, 0.68859256, 0.558599,0.56255131,0.61822723,\
                                        0.58297106, 0.67802997, 0.56411052,0.70597755]
                        if name == 'u603': 
                            meanpear = [0.60075925, 0.62926192, 0.7419749,0.69909184, 0.63475069, 0.52538688,0.61374631, 0.61529876, 0.53526537,\
                                        0.61416945, 0.515072,0.74873393]
                        self.pear0 = np.append(self.pear0,meanpear[0])
                        self.pear1 = np.append(self.pear1,meanpear[1]) 
                        self.pear2 = np.append(self.pear2,meanpear[2]) 
                        self.pear3 = np.append(self.pear3,meanpear[3]) 
                        self.pear4 = np.append(self.pear4,meanpear[4]) 
                        self.pear5 = np.append(self.pear5,meanpear[5]) 
                        self.pear6 = np.append(self.pear6,meanpear[6])  
                        self.pear7 = np.append(self.pear7,meanpear[7]) 
                        self.pear8 = np.append(self.pear8,meanpear[8]) 
                        self.pear9 = np.append(self.pear9,meanpear[9]) 
                        self.pear10 = np.append(self.pear10,meanpear[10]) 
                        self.pear11 = np.append(self.pear11,meanpear[11]) 
                        if name == 'u601' or name == 'u602':
                            self.pear12 = np.append(self.pear12,meanpear[12]) 
                        if name == 'u601':
                            self.pear13 = np.append(self.pear13,meanpear[13])  

                self.pears = [self.pear0,self.pear1,self.pear2,self.pear3,self.pear4,self.pear5,self.pear6,\
                             self.pear7,self.pear8,self.pear9,self.pear10,self.pear11,self.pear12,self.pear13]
                self.pearsonr = np.empty([0],dtype=float)  
      

        # - - - - - - - - - - - - - - - - - - - - 
        # CONTINUE RUN DEFINITION  
        all_cores = np.unique(thtr.core_ids) 
        if core_list is None:
            core_list = all_cores 
        # - - - - - - - - - - - - - - - - - - - - 
        # CORELOOP
        for nc,core_id in enumerate(core_list):  
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=False) 
            tmap=rainbow_map(ms.ntimes) 
 
            # FIELDS
            rho_extents=davetools.extents()
            magfield_extents = davetools.extents()  

            density = thtr.c([core_id],'density') 
            magfield = thtr.c([core_id],'magnetic_field_strength') 
            cv = thtr.c([core_id],'cell_volume')  
 
            rho_extents(density)
            magfield_extents(magfield) 
           
            asort =  np.argsort(thtr.times)
            if (asort != sorted(asort)).any():
                print("Warning: times not sorted.") 

            # FIG 1.
            if typeplot == 'scatter_plot':
                #fig, ax1=plt.subplots(1,1) 
                if core_id == core_val:
                    fig, (ax1, ax2) = plt.subplots(1, 2)
                    fig.subplots_adjust(wspace=0, hspace=0)
            # - - - - - - - - - - - - - - - - - - - - 
            # TIMELOOP: 
            densitymasked = np.empty([0],dtype=bool)
            magfieldmasked = np.empty([0],dtype=bool)

            for n_count,n_time in enumerate(asort):       
                mask = ms.compute_unique_mask(core_id, dx=1./2048,frame=n_time)  
                c=tmap(n_count,mask.sum()) 
                
                if timing == 'per_frame':    
                    XX = np.log10(density[mask,n_time])
                    X = np.log10(density[mask,n_time]).flatten()  # [mask,n_time], OR [:,n_time][mask]
                    X2 = np.linspace(X.min(),X.max(),num=len(X))
                    XX2 = 10 ** X2 

                    YY = np.log10(magfield[mask,n_time])
                    Y = np.log10(magfield[mask,n_time]).flatten() 
                    C = np.log(cv[mask,n_time]).flatten()

                if timing == 'all_time':
                    densityflat = np.log10(density[mask,n_time]).flatten()
                    magfieldflat = np.log10(magfield[mask,n_time]).flatten()
                
                    densitymasked = np.append(densitymasked,densityflat)
                    magfieldmasked = np.append(magfieldmasked, magfieldflat)

                if timing == 'all_time' and n_time == asort[-1]:
                    X = densitymasked
                    X2 = np.linspace(X.min(),X.max(),num=len(X)) 
                    XX2 = 10 ** X2 
                    Y = magfieldmasked

                elif timing == 'all_time': 
                    XX2 = 0
                    YY = 0
                    YY_mean = 0
                               
                # - - - - - - - - - - - - - - - - - - - - 
                # POLY-FITTING:
                if timing == 'per_frame' or n_time == asort[-1]:
                    pfit = np.polyfit(X,Y,1) 
                    other = np.poly1d(pfit)
                    beta = pfit[0]
                    B_o = pfit[1]  
                    YY = 10 ** (pfit[0]*X2 + pfit[1]) 
                    self.betarr = np.append(self.betarr,beta)                         
                    
                    # - - - - - - - - - - - - - - - - - - - - 
                    #PEARSON FOR ALL TIME
                    if typeplot == 'vio_plot':
                        allxs = np.std(X)  
                        allys = np.std(Y) 
                        if allxs != 0 and allys != 0:
                            pearX,pearY = scipy.stats.pearsonr(X,Y)
                            self.pearsonr = np.append(self.pearsonr,pearX) 
                        else:
                            print("A zero encountered!!",allxs,allys)
 
                    # plot the alpha of the mean...manually for now 
                    # these are mask-ed
                    if name == 'u601':
                        YY_mean = 10 ** ((0.204)*X2 + pfit[1])
                    if name == 'u602':
                        YY_mean = 10 ** ((0.313)*X2 + pfit[1])
                    if name == 'u603':
                        YY_mean = 10 ** ((0.373)*X2 + pfit[1])

                # FIGs 1, 3, 5, 6 
                # CALL PEARSON R
                if typeplot == 'vio_plot':   
                    pearsonR(X,Y)
                    scatterplots()  
                # CALL SCATTER:
                if typeplot == 'frame_scatters':
                    scatterplots(fig,ax1,ax2,ax3,ax4,ax5,ax7,ax8,lplots,XX2,YY)  
                if typeplot == 'scatter_plot':  
                    scatterplots(fig,ax1,ax2,xx2=XX2,yy=YY,yym=YY_mean)  
                if typeplot == 'box_plot':
                    scatterplots()
           
            plt.close(fig)  # at this position, †his closes the plot of each core for all time
        return ax2,fig


    def twopanels(self,num,name,ax2,fig,coreval):
        print('inside two panels')
        # this panel needs to run ALL CORES
        the_min = self.betarr.min()
        the_max = self.betarr.max()
        the_bins = math.isqrt(len(self.betarr))
        the_FTA = np.mean(self.betarr)
        #atf = [0.175,0.270,0.316]
        #the_ATF = atf[num]  
        the_bins = np.linspace(the_min,the_max,num=the_bins)  #'original': 64
        ax2.hist(self.betarr, bins=the_bins, density=True, histtype='step', color='k')
        ax2.axvline(the_FTA, color='grey', linestyle='dashed')
        #ax2.axvline(the_ATF, color='m')

        ax2.yaxis.tick_right()
        ax2.set_xlabel(r'$\alpha$') 
        ax2.set_ylabel('%') 
        y_vals = ax2.get_yticks()
        ax2.set_yticklabels(['{:.3f}'.format(x/len(self.betarr)) for x in y_vals])

        outname = 'brhotff_histos_c%04d_%s'%(coreval,name)  
        fig.savefig(outname) 
        print("saved "+outname)

            
    def boxes(self,num,name,tff_lab): 
        Bears = {}
        index = []
        data = [] 
 
        for i,tff in enumerate(tff_lab): 
            Bears[tff_lab[i]] = self.bears[i]
        for j, (key, val) in enumerate(Bears.items()):
            index.append(key)
            data.append(val)
        fig, ax1=plt.subplots(1,1)    

        meanpointprops = dict(marker='.', markeredgecolor='black',markerfacecolor='blue') 
        bparts = ax1.boxplot(data, meanprops=meanpointprops, showmeans=True, meanline=False, showfliers=False)
        for i,val in enumerate(tff_lab):
            y = data[i] 
            x = np.random.normal(1+i, 0.04, size=len(y))
            ax1.plot(x, y, 'b.', alpha=0.2)
            if 1:
                y2 = data[i][-1]
                x2 = x[-1] 
                ax1.plot(x2,y2,'r.',alpha=1)  
                                          #CAREFUL I have now changed the mean ever so slightly by adding the meanbears :\
 
        a = 0.4
        if num == 0: 
            meanbear = [0.207432,0.201813,0.234106,0.230126,0.160165,0.074837,0.107986,0.079870,0.030422,0.043188,0.101638,0.176109,0.353048,0.456282]    
            ax1.plot([1,2,3,4,5,6,7,8,9,10,11,12,13,14],meanbear,c='r',alpha=a)
            ax1.plot([1,14],[0,0],'--',c=[0.5]*4) 
            ax1.plot([1,14],[0.667,0.667],'--',c=[0.5]*4)
            ax1.plot([1,14],[0.5,0.5],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.07,0.15,0.23,0.30,0.38,0.45,0.53,0.61,0.68,0.75,0.82,0.89,0.95] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1))
            title = r'$\beta = 0.2$'
 
        if num == 1: 
            meanbear = [0.273494,0.265940,0.329641,0.385815,0.344329,0.284348,0.203154,0.173078,0.160007,0.166019,0.240440,0.277389,0.479462]
            ax1.plot([1,2,3,4,5,6,7,8,9,10,11,12,13],meanbear,c='r',alpha=a)
            ax1.plot([1,13],[0,0],'--',c=[0.5]*4) 
            ax1.plot([1,13],[0.667,0.667],'--',c=[0.5]*4)
            ax1.plot([1,13],[0.5,0.5],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.08,0.16,0.25,0.33,0.41,0.50,0.58,0.66,0.74,0.82,0.90,0.97] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1))
            title = r'$\beta = 2.0$'

        if num == 2: 
            meanbear = [0.279555,0.336840,0.355424,0.325673,0.268460,0.232444,0.258084,0.252883,0.206200,0.239367,0.283474,0.617871]
            ax1.plot([1,2,3,4,5,6,7,8,9,10,11,12],meanbear,c='r',alpha=a)
            ax1.plot([1,12],[0,0],'--',c=[0.5]*4) 
            ax1.plot([1,12],[0.667,0.667],'--',c=[0.5]*4)
            ax1.plot([1,12],[0.5,0.5],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.08,0.17,0.26,0.35,0.44,0.52,0.60,0.69,0.77,0.86,0.91] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1)) 
            title = r'$\beta = 20$'

        ax1.set_xticklabels(tff_mod)       
        ylim = -0.5,1.0
        BRho_tool.labelled(ax1,xlabel=r'$t/t_{\rm{ff}}$',ylabel=r'$\alpha$',ylim=ylim) 
        fig.savefig('BearsBoxplot_cf%s'%name)
        print("BOX SAVED")
        plt.close(fig)
      

    def violins(self,num,name,tff_lab):
    # NOTE: for PearsonR, nan values are not accepted 
        Pears = {}
        index = []
        data = []  

        for i,tff in enumerate(tff_lab):
            Pears[tff_lab[i]] = self.pears[i] 
        for j, (key, val) in enumerate(Pears.items()):
            index.append(key)
            data.append(val)

        fig, ax1=plt.subplots(1,1)             

        meanpointprops = dict(marker='.', markeredgecolor='black',markerfacecolor='blue') 
        bparts = ax1.boxplot(data,meanprops=meanpointprops,showmeans=True,meanline=False,showfliers=False)
        for i,val in enumerate(tff_lab):
            y = data[i] 
            x = np.random.normal(1+i, 0.04, size=len(y))
            ax1.plot(x, y, 'b.', alpha=0.2)
            if 1:
                y2 = data[i][-1]
                x2 = x[-1] 
                ax1.plot(x2,y2,'r.',alpha=1)
                
        #To compare violinplots with a zero; affects the position of the plots
        if num == 0:
            ax1.plot([1,14],[0,0],c=[0.5]*4)  
            ax1.plot([1,14],[1,1],'--',c=[0.5]*4)
            ax1.plot([1,14],[-1,-1],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.07,0.15,0.23,0.30,0.38,0.45,0.53,0.61,0.68,0.75,0.82,0.89,0.95] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1))
            
        if num == 1:
            ax1.plot([1,13],[0,0],c=[0.5]*4) 
            ax1.plot([1,13],[1,1],'--',c=[0.5]*4)
            ax1.plot([1,13],[-1,-1],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.08,0.16,0.25,0.33,0.41,0.50,0.58,0.66,0.74,0.82,0.90,0.97] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1))

        if num == 2:
            ax1.plot([1,12],[0,0],c=[0.5]*4) 
            ax1.plot([1,12],[1,1],'--',c=[0.5]*4)
            ax1.plot([1,12],[-1,-1],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.08,0.17,0.26,0.35,0.44,0.52,0.60,0.69,0.77,0.86,0.91] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1)) 
       
        ax1.set_xticklabels(tff_mod)      
        ylim = -1.15,1.15 
        BRho_tool.labelled(ax1,xlabel=r'$t/t_{\rm{ff}}$',ylabel=r'$R$',ylim=ylim)
         
        fig.savefig('PearsViolin_%s'%name)
        print("SAVED")
        plt.close(fig)

# - - - - - - - - - - - - - - - - - - - - 
# "MAIN"
# THREE SIMS AT ONCE
import three_loopers_six as TL6
if 'clobber' not in dir():
    clobber=True

if 'BRho_tool1' not in dir() or clobber: 
    BRho_tool1=BRho_tool(TL6.loops['u601'])
    simname1 = 'u601' 

if 'BRho_tool2' not in dir() or clobber: 
    BRho_tool2=BRho_tool(TL6.loops['u602']) 
    simname2 = 'u602'
  
if 'BRho_tool3' not in dir() or clobber: 
    BRho_tool3=BRho_tool(TL6.loops['u603'])  
    simname3 = 'u603'
simnames = [simname1,simname2,simname3]
print("GREETINGS")

'''
PREP YOUR FIGURE:
FIGURE 1: type: 'scatter_plot' ;time: 'all_time'

FIGURE 3: type: 'frame_scatters' ;time: 'per_frame'
    
FIGURE 5: type: 'box_plot' ;time: 'per_frame'

FIGURE 6: type: 'vio_plot' ;time: 'per_frame'
'''

for nt,tool in enumerate([BRho_tool1,BRho_tool2,BRho_tool3]):
    # TYPE: 
    which_plot = 'scatter_plot' 
    # TIME:
    which_time = 'all_time'

    # GLOBAL TFF 
    G = 1620/(4*np.pi)
    rho_mean = 1
    t_ff = np.sqrt(3*np.pi/(32*G*rho_mean)) 

    # TFF PERCENTAGE & DESIRED FRAMES, DESIRED CORES
    # fix "core_list" respectively! 
    # corenum: representative core per simulation
    thtr = tool.this_looper.tr
    if nt == 0: 
        tff_p = thtr.times/t_ff 
        frames = [1,3,6,8,10,12]   
        core_list = [27]
        corenum = 27
    if nt == 1:
        tff_p = thtr.times/t_ff 
        frames = [1,3,5,7,9,11] 
        core_list = [32]
        corenum = 32
    if nt == 2:
        tff_p = thtr.times/t_ff 
        frames = [1,3,5,7,9,10] 
        core_list = [98]
        corenum = 98
    
    simframes = set(frames) 
    tff_labels = ['%.2f'%s for s in tff_p]
    
    if which_plot == 'frame_scatters': 
        fig = plt.figure() 
        fig.text(0.365,0.03,r'$\rho/\rho_{o}$')

        ax1 = plt.subplot(331)
        ax2 = plt.subplot(332)
        ax4 = plt.subplot(334)
        ax5 = plt.subplot(335)
        ax7 = plt.subplot(337)
        ax8 = plt.subplot(338)

        ax3 = plt.subplot(133)  #`ax6, ax9` 
        fig.subplots_adjust(wspace=0, hspace=0)
        
        if nt == 0:
            lplots = [0,ax1,0,ax2,0,0,ax4,0,ax5,0,ax7,0,ax8,0] 
        if nt == 1:
            lplots = [0,ax1,0,ax2,0,ax4,0,ax5,0,ax7,0,ax8,0] 
        if nt == 2:
            lplots = [0,ax1,0,ax2,0,ax4,0,ax5,0,ax7,ax8,0] 

        tool.run(simnames[nt],thtr,which_time,which_plot,tff_labels,
             fig,ax1,ax2,ax3,ax4,ax5,ax7,ax8,
             lplots,core_list,core_val=None,simframes=simframes) 

    else:  
        # assign 'None' respectively 
        print('about to run!')
        axnum,figu = tool.run(simnames[nt],thtr,which_time,which_plot,tff_labels,
             fig=None,ax1=None,ax2= None,ax3= None,ax4= None,ax5= None,ax7= None,ax8= None,
             lplots= None,core_list=None,core_val=corenum,simframes= None) 

    if which_plot == 'scatter_plot' and which_time == 'all_time': 
        # TWO PANEL SCATTER PLOTS; Fig 1.
        tool.twopanels(nt, simnames[nt],axnum,figu,corenum)
    if which_plot == 'box_plot' and which_time == 'per_frame': 
        # CALL BOXPLOTS; Fig 5. 
        tool.boxes(nt,simnames[nt],tff_labels)
    if which_plot == 'vio_plot' and which_time == 'per_frame': 
        # CALL VIOLINPLOTS; Fig 6.              
        tool.violins(nt,simnames[nt],tff_labels)





