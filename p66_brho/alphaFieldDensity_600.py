
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


plt.close('all')


# - - - - - - - - - - - - - - - - - - - - 
# THE CLASS
class BRho_tool():
    def __init__(self,this_looper): 
        self.this_looper = this_looper
        self.betarr = np.empty([0],dtype=float)
        self.pearsonr = np.empty([0],dtype=float)
        self.corestdzero = np.empty([0],dtype=float)

        # need to add 'self' to the following in the body of code: 
        self.time_stamp = np.empty([0],dtype=float) 
        self.pearson_lmr = np.empty([0],dtype=float)
        self.pearson_mr = np.empty([0],dtype=float)
        self.spearmanr = np.empty([0],dtype=float)
   
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

    # CONFIRM THIS WORKS 
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

        print('core list',core_list)
        print('simframes are ',simframes)     
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
                #print("core id and sim ",core_id,name)  #EXP
                #self.corestdzero = np.append(self.corestdzero,core_id)  #EXP
                pearX = 0  # try a break from the current loop...which would only be of the time, not core... 
            self.pearsonr = np.append(self.pearsonr,pearX)  
                

        def scatterplots(fig=None,ax1=None,ax2=None,ax3=None,ax4=None,ax5=None,ax7=None,ax8=None,
                         lplots=None,xx2=None,yy=None,yym=None):
           
            if typeplot == 'scatter_plot' and core_id == core_val:
                if n_time <= asort[-2]:
                    #print('n_time',n_time)  #TIME CHECK
                    ax1.scatter(density[mask,n_time],magfield[mask,n_time],c=c,label=thtr.times[n_time],s=0.1)  #edit thtr.times            

                xlims = 10e-3,10e6      
                ylims = 10e-2,10e3      
                if n_time == asort[-2]:  #this way THIS DOES NOT SAVE THE FINAL erronous frame  
                    #print('tfflabs for printing',tfflabs[n_time])  #TIME CHECK
                    ax1.plot(xx2,yy,c='k',linewidth=1.0) #c=c[0] for colors
                    ax1.plot(xx2,yym,'--',c='k',linewidth=1.0) #c=c[0] for colors
                    BRho_tool.labelled(ax1,xscale='log',yscale='log',xlabel=r'$\rho/\rho_{o}$',ylabel=r'$\mid B \mid (\mu G)$',
                                       xlim=xlims, ylim=ylims)#,title=r'$\alpha_{b\rho} = %.3f$'%beta)
                                       #xlim=rho_extents.minmax, ylim=magfield_extents.minmax,title=r'$\alpha_{b\rho} = %.3f$'%beta) 

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
                        ax3.set_title(r"$\beta = 0.2, \alpha_{cf}$")
                    if name == 'u602':
                        the_range = np.arange(0.075,1,0.075) 
                        ax3.scatter(the_range[:-1],self.betarr[:-1],c=c2)  
                        ax3.plot(the_range[:-1],self.betarr[:-1],c='k',linewidth=1.0)  
                        ax3.set_title(r"$\beta = 2.0, \alpha_{cf}$")
                    if name == 'u603':
                        the_range = np.arange(0.075,0.9,0.075) 
                        ax3.scatter(the_range[:-1],self.betarr[:-1],c=c2)  
                        ax3.plot(the_range[:-1],self.betarr[:-1],c='k',linewidth=1.0)   
                        ax3.set_title(r"$\beta = 20, \alpha_{cf}$")
                    ax3.plot([0,1],[0,0],c=[0.5]*4) 
                    ax3.plot([0,1],[0.4,0.4],'--',c=[0.5]*4)
                 
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
                if name == 'u601':
                    self.bear12 = np.append(self.bear12,self.betarr[12])    
                    self.bear13 = np.append(self.bear13,self.betarr[13]) 

                if 1:
                    if core_id == core_list[-1]:
                        if name == 'u601':
                            meanbear = [0.207432,0.201813,0.234106,0.230126,0.160165,0.074837,0.107986,0.079870,0.030422,0.043188,0.101638,0.176109,0.353048,0.456282]    
                            fitbear = [0.037122767966842135,0.06130065264472109,0.10876453322815462,0.12037348290995548,0.1256151386016793,0.12584476500389483,\
                                       0.12960541113191903,0.14452750478517815,0.16035081563879608,0.17350352881649694,0.20320963091004954,0.22354626336366118,\
                                       0.27401985393584366,0.3323051356106275]
                        if name == 'u602':
                            meanbear = [0.273494,0.265940,0.329641,0.385815,0.344329,0.284348,0.203154,0.173078,0.160007,0.166019,0.240440,0.277389,0.479462]
                            fitbear = [0.15059829575793074,0.19999384737531586,0.24022574462120863,0.27614753632431754,0.27614753632431754,0.29332289603778755,\
                                       0.3097313041491559,0.30801461152950443,0.31782695369769487,0.3358170167076137,0.3148400500596815,0.28845754861915246,\
                                       0.41477394063271383]
                        if name == 'u603': 
                            meanbear = [0.279555,0.336840,0.355424,0.325673,0.268460,0.232444,0.258084,0.252883,0.206200,0.239367,0.283474,0.617871]
                            fitbear = [0.22829863399090444,0.3185318136356117,0.37526692580639975,0.40111503807258014,0.40640704584839654,0.38681096100128876,\
                                       0.38332517520463627,0.39490446738439916,0.371243675696434,0.3321847915857663,0.32833900753678574,0.3893047090450924]
                        self.bear0 = np.append(self.bear0,(meanbear[0],fitbear[0]))
                        self.bear1 = np.append(self.bear1,(meanbear[1],fitbear[1]))
                        self.bear2 = np.append(self.bear2,(meanbear[2],fitbear[2])) 
                        self.bear3 = np.append(self.bear3,(meanbear[3],fitbear[3])) 
                        self.bear4 = np.append(self.bear4,(meanbear[4],fitbear[4])) 
                        self.bear5 = np.append(self.bear5,(meanbear[5],fitbear[5])) 
                        self.bear6 = np.append(self.bear6,(meanbear[6],fitbear[6]))  
                        self.bear7 = np.append(self.bear7,(meanbear[7],fitbear[7])) 
                        self.bear8 = np.append(self.bear8,(meanbear[8],fitbear[8])) 
                        self.bear9 = np.append(self.bear9,(meanbear[9],fitbear[9]))
                        self.bear10 = np.append(self.bear10,(meanbear[10],fitbear[10])) 
                        self.bear11 = np.append(self.bear11,(meanbear[11],fitbear[11]))
                        if name == 'u601':
                            self.bear12 = np.append(self.bear12,(meanbear[12],fitbear[12]))    
                            self.bear13 = np.append(self.bear13,(meanbear[13],fitbear[13])) 

                self.bears = [self.bear0,self.bear1,self.bear2,self.bear3,self.bear4,self.bear5,self.bear6,\
                              self.bear7,self.bear8,self.bear9,self.bear10,self.bear11,self.bear12,self.bear13]

                self.betarr = np.empty([0],dtype=float)  


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
                self.pear11 = np.append(self.pear11,self.pearsonr[11])  #final u603 pear contains nans 
                if name == 'u601':
                    self.pear12 = np.append(self.pear12,self.pearsonr[12])  #final u602 pear, but want to discard final frame   
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
                        self.pear11 = np.append(self.pear11,meanpear[11])  #final u603 pear contains nans 
                        if name == 'u601':
                            self.pear12 = np.append(self.pear12,meanpear[12])  #final u602 pear, but want to discard final frame   
                            self.pear13 = np.append(self.pear13,meanpear[13])  

                self.pears = [self.pear0,self.pear1,self.pear2,self.pear3,self.pear4,self.pear5,self.pear6,\
                             self.pear7,self.pear8,self.pear9,self.pear10,self.pear11,self.pear12,self.pear13]
                self.pearsonr = np.empty([0],dtype=float)  
       
        # - - - - - - - - - - - - - - - - - - - - 
        # CONTINUE RUN DEFINITION  
        all_cores = np.unique(thtr.core_ids) 
        if core_list is None:
            core_list = all_cores
 
        #rm = rainbow_map(len(core_list))  # or all_cores, needed??  
        rho_extents=davetools.extents()
        magfield_extents = davetools.extents()  


        # - - - - - - - - - - - - - - - - - - - - 
        # CORELOOP
        for nc,core_id in enumerate(core_list):  
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=False) 
            tmap=rainbow_map(ms.ntimes) 
 
            # FIELDS
            density = thtr.c([core_id],'density') 
            magfield = thtr.c([core_id],'magnetic_field_strength')
            magden = magfield/density
            cv = thtr.c([core_id],'cell_volume')  
 
            if name != 'u601': 
                B_x = thtr.c([core_id],'magnetic_field_x')  
                B_y = thtr.c([core_id],'magnetic_field_y')
                B_z = thtr.c([core_id],'magnetic_field_z')
                v_x = thtr.c([core_id],'velocity_x')
                v_y = thtr.c([core_id],'velocity_y')
                v_z = thtr.c([core_id],'velocity_z') 

            rho_extents(density)
            magfield_extents(magfield) 
           
            asort =  np.argsort(thtr.times)
            if (asort != sorted(asort)).any():
                print("Warning: times not sorted.") 

            # SETUP FIGURES FOR EACH CORE HERE for all time: 
            if typeplot == 'scatter_plot':
                #fig, ax1=plt.subplots(1,1) 
                if core_id == core_val:
                    fig, (ax1, ax2) = plt.subplots(1, 2)
                    fig.subplots_adjust(wspace=0, hspace=0)
            if typeplot == 'rms_plot':
                fig, ax1=plt.subplots(1,1)
           

            # - - - - - - - - - - - - - - - - - - - - 
            # TIMELOOP: 
            densitymasked = np.empty([0],dtype=bool)
            magfieldmasked = np.empty([0],dtype=bool)
            for n_count,n_time in enumerate(asort):       
                mask = ms.compute_unique_mask(core_id, dx=1./2048,frame=n_time)  
                #print('mask shape',mask.shape)
                c=tmap(n_count,mask.sum())  #Used to be: ms.nparticles 
                
                if timing == 'per_frame':    
                    #print('n_time',n_time)  #TIME CHECK
                    XX = np.log10(density[mask,n_time])
                    X = np.log10(density[mask,n_time]).flatten()  # [mask,n_time], OR [:,n_time][mask]
                    X2 = np.linspace(X.min(),X.max(),num=len(X)) # OR PER FRAME, MULTIPLE PANEL PLOTS 
                    XX2 = 10 ** X2 

                    YY = np.log10(magfield[mask,n_time])
                    Y = np.log10(magfield[mask,n_time]).flatten() 

                    Z = np.log10(magden[mask,n_time]).flatten()
                    C = np.log(cv[mask,n_time]).flatten()
                    WW = YY/XX
                    W = WW.flatten() 

                if timing == 'all_time':
                    #time=thtr.times[n_time]
                    #print('tfflabs to append masks',tfflabs[n_time])  #TIME CHECK

                    densityflat = np.log10(density[mask,n_time]).flatten()
                    magfieldflat = np.log10(magfield[mask,n_time]).flatten()
                
                if 0:  #FOR ALL TIME ALPHA
                    densitymasked = np.append(densitymasked,densityflat)
                    magfieldmasked = np.append(magfieldmasked, magfieldflat)

                if 0:  #FOR CRUTCHER'S HD ALPHA
                    if densityflat.max() > 2.602059991:
                        densitymasked = np.append(densitymasked,densityflat)
                        magfieldmasked = np.append(magfieldmasked, magfieldflat)

                if timing =='all_time' and n_time == asort[-2]:
                    #print('tfflabs masked for polyfit',tfflabs[n_time])  #TIME CHECK
                    X = densitymasked
                    X2 = np.linspace(X.min(),X.max(),num=len(X))  # FOR ALL TIME, ONE PANEL PLOTS 
                    XX2 = 10 ** X2 
                    Y = magfieldmasked

                elif timing == 'all_time': 
                    #print('tfflabs wait while being masked',tfflabs[n_time])  #TIME CHECK
                    XX2 = 0
                    YY = 0
                    YY_mean = 0
                               
                # POLY-FITTING:
                if timing == 'per_frame' or n_time == asort[-2]:
                    #print('tfflabs the polyfit',tfflabs[n_time])  #TIME CHECK
                    pfit = np.polyfit(X,Y,1) 
                    other = np.poly1d(pfit)
                    beta = pfit[0]
                    B_o = pfit[1]  
                    YY = 10 ** (pfit[0]*X2 + pfit[1]) 

                    self.betarr = np.append(self.betarr,beta)                         
                    # THE NEGATIVE OUTLIERS
                    #if beta < 0: 
                    #    self.time_stamp = np.append(self.time_stamp,n_time)

                    #PEARSON FOR ALL TIME
                    if 0:
                        allxs = np.std(X)  
                        allys = np.std(Y) 
                        if allxs != 0 and allys != 0:
                            pearX,pearY = scipy.stats.pearsonr(X,Y)
                            self.pearsonr = np.append(self.pearsonr,pearX) 
                        else:
                            print("A zero encountered!!",allxs,allys)
 
                    # plot the alpha of the mean;  #may need to be edited 
                    # these have mask & not the last frame incorporated
                    if name == 'u601':
                        YY_mean = 10 ** ((0.196)*X2 + pfit[1])
                    if name == 'u602':
                        YY_mean = 10 ** ((0.314)*X2 + pfit[1])
                    if name == 'u603':
                        YY_mean = 10 ** ((0.383)*X2 + pfit[1])

 
                # FOR HISTOGRAMS OF DESIRED VALUES
                if 0:
                    x = density[mask,n_time]
                    y = magfield[mask,n_time]
                    self.bears[n_count] = np.append(self.bears[n_count],Z) 
                    self.pears[n_count] = np.append(self.pears[n_count],W) 
                    self.dears[n_count] = np.append(self.dears[n_count],C) 
              

                # CALL PEARSON R
                if typeplot == 'vio_plot':   
                    pearsonR(X,Y)
                    scatterplots()  
                # CALL SCATTER: ...unfortunately "or" behaves strangely
                if typeplot == 'frame_scatters':
                    scatterplots(fig,ax1,ax2,ax3,ax4,ax5,ax7,ax8,lplots,XX2,YY)  
                if typeplot == 'scatter_plot':  
                    scatterplots(fig,ax1,ax2,xx2=XX2,yy=YY,yym=YY_mean)  
                if typeplot == 'box_plot':
                    scatterplots()


                # RMS EXPERIMENTS
                if typeplot == 'rms_plot':  
                    if name != 'u601':  
                        B_avg = (magfield[mask,n_time] * cv[mask,n_time]).sum()/cv[mask,n_time].sum()  

                        Bx_avg = (B_x[mask,n_time] * cv[mask,n_time]).sum()/cv[mask,n_time].sum() 
                        By_avg = (B_y[mask,n_time] * cv[mask,n_time]).sum()/cv[mask,n_time].sum()
                        Bz_avg = (B_z[mask,n_time] * cv[mask,n_time]).sum()/cv[mask,n_time].sum()
                        vx_avg = (v_x[mask,n_time] * cv[mask,n_time]).sum()/cv[mask,n_time].sum()
                        vy_avg = (v_y[mask,n_time] * cv[mask,n_time]).sum()/cv[mask,n_time].sum()
                        vz_avg = (v_z[mask,n_time] * cv[mask,n_time]).sum()/cv[mask,n_time].sum()

                        Bx_sq = (B_x[mask,n_time] - Bx_avg)**2 
                        By_sq = (B_y[mask,n_time] - By_avg)**2
                        Bz_sq = (B_z[mask,n_time] - Bz_avg)**2
                        vx_sq = (v_x[mask,n_time] - Bx_avg)**2
                        vy_sq = (v_y[mask,n_time] - By_avg)**2
                        vz_sq = (v_z[mask,n_time] - Bz_avg)**2
                    
                        # SPREAD OF B FIELD FROM THE MEAN, also the standard dev...these should give same-sh as np.std()?  
                        brms = np.sqrt(np.mean(Bx_sq) + np.mean(By_sq) + np.mean(Bz_sq))  # DO MEAN WITH CV !
                        brms_py = np.std(magfield)  #COMPARE for fun, but for the vector field
                        #print('BRMS',brms)
                        #print('BRMS_py',brms_py)
                        vrms = np.sqrt(np.mean(vx_sq) + np.mean(vy_sq) + np.mean(vz_sq))  # DO MEAN WITH CV!! 
                        #vrms_py = np.std()   #COMPARE for fun, VELOCITY mag.
                        vrmsq_bmean = (vrms**2)/B_avg  #B_avg for all box or region; try diff ways
                        #print('VRMS',vrms)
                        #print('VRMS_py',vrms_py)

                        self.brms = np.append(self.brms,brms) 
                        self.vrms = np.append(self.vrms,vrms)
                        self.vrmsq_bmean = np.append(self.vrmsq_bmean,vrmsq_bmean)
                        
                        if n_time == asort[-1]:
                            print('vrmsq_bmean',self.vrmsq_bmean)
                            #breakpoint()
                            print('inside!')
                            tmap2 = rainbow_map(len(self.brms)) 
                            c2 = [tmap2(n) for n in range(len(self.brms))]  

                            if name == 'u602':
                                the_range = np.arange(0.075,1,0.075)  
                                #ax1.scatter(the_range,self.brms,c='b')  
                                #ax1.scatter(the_range,self.vrms,c='g')  
                                #ax1.scatter(self.vrms, self.brms,c=c2) 
                                ax1.scatter(self.vrmsq_bmean, self.brms,c=c2) 
                                #ax1.set_xscale('log')
                                #ax1.set_yscale('log')
                                #ax1.set_title(r"$\beta = 2.0$, B_rms: blue, V_rms: green, core:%04d"%core_id)
                                ax1.set_title(r"$\beta = 2.0$, B_rms vs V_rms vs tff, core:%04d"%core_id)
                                
                                outname = 'Brms_vrmsqBm_tff_c%04d_%s'%(core_id,name)  
                                plt.savefig(outname)  # CAREFUL WITH FIG VS PLT
                                print("saved "+outname)
                                #plt.close(fig)  #CAREFUL, this doesn't allow to do another desired core :(   
                                
                            if name == 'u603':
                                the_range = np.arange(0.075,0.9,0.075) 
                                #ax1.scatter(the_range,self.brms,c='b')  
                                #ax1.scatter(the_range,self.vrms,c='g')  
                                #ax1.scatter(self.vrms, self.brms,c=c2)  
                                ax1.scatter(self.vrmsq_bmean, self.brms,c=c2) 
                                #ax1.set_xscale('log')
                                #ax1.set_yscale('log')
                                #ax1.set_title(r"$\beta = 20.$, B_rms: blue, V_rms: green, core:%04d"%core_id)
                                ax1.set_title(r"$\beta = 20.$, B_rms vs V_rms vs tff, core:%04d"%core_id)

                                outname = 'Brms_vrmsqBm_tff_c%04d_%s'%(core_id,name)  
                                plt.savefig(outname)  # CAREFUL WITH FIG VS PLT
                                print("saved "+outname)
                                #plt.close(fig)  #CAREFUL, this doesn't allow to do another desired core :(   
                            
                            self.vrmsq_bmean = np.empty([0],dtype=float)                               
                            self.brms = np.empty([0],dtype=float)                               
                            self.vrms = np.empty([0],dtype=float)                               

           
            plt.close(fig)  # at this position, †his closes the plot of each core for all time
            if 0:  # write alpha means of all cores for all time at once
                if core_id == core_list[-1]:
                    alphaFile = open("alphaRecords.txt",'a')
                    # the following 3 prints should match
                    print("core_list_len",len(core_list))
                    print("lenAlphas ",len(self.betarr))   
                    print("lenPears ",len(self.pearsonr))   
                    meanAlpha = np.mean(self.betarr)
                    meanPear = np.mean(self.pearsonr)
                    print('meanAlpha ',meanAlpha)
                    print('meanPear',meanPear)
                    alphaFile.write("Sim %s meanAlpha %f \n"%(name,meanAlpha))
                    alphaFile.write("Sim %s meanPear %f \n"%(name,meanPear))
                    alphaFile.close()
        return ax2,fig


    def twopanels(self,num,name,ax2,fig,coreval):
        # this panel needs to run ALL CORES
        ax2.hist(self.betarr, 64, density=True, histtype='step', color='k')
        ax2.yaxis.tick_right()
        ax2.set_xlabel(r'$\alpha$') 
        ax2.set_ylabel('%') 
        y_vals = ax2.get_yticks()
        ax2.set_yticklabels(['{:.3f}'.format(x/len(self.betarr)) for x in y_vals])

        outname = 'brhotff_histos_c%04d_plmean_%s'%(coreval,name)  
        fig.savefig(outname)  # CAREFUL WITH FIG VS PLT
        print("saved "+outname)

    def histograms(self,num,figs=None,plts=None):
    # HISTOGRAM FOR ALL CORES FOR ALL TIME, now per frame.

        #BETA ENTRIES, AVERAGE & STD 
        betavg = np.mean(self.betarr)
        betastd = np.std(self.betarr)
        entries = len(self.betarr) 
     
        # EDIT: I should do RGB :p
        if num == 0:
            color = 'r'
        if num == 1:
            color = 'b'
        if num == 2:
            color = 'g'

        #plts.hist(self.betarr, 50, density=False, histtype='step', color=color)
        for i in range(len(self.bears)):

            if 1:
                # HISTOS PER FRAME: EDITTTTTT 
                the_bins = np.linspace(-10,10)      #need to find min and max of bears and pearsi, or same as mean
                the_weight = self.dears[i]          #need to make self.dears[i];MADE, TRY; try weights and no weights,correct y axis accordingly

                # THE LN(Y/X)
                the_lnArray, xbins = np.histogram(self.bears[i],bins=the_bins,weights=the_weight,density=True)
                bin_lncenters = 0.5*(xbins[1:]+xbins[:-1])
                the_lnX= bin_lncenters
                the_lnY= the_lnArray
                # THE LNY/LNX
                the_lnlnArray, xxbins = np.histogram(self.pears[i],bins=the_bins,weights=the_weight,density=True)
                bin_lnlncenters = 0.5*(xxbins[1:]+xxbins[:-1])
                the_lnlnX= bin_lnlncenters
                the_lnlnY= the_lnlnArray

                plts.plot(the_lnX,the_lnY,c=color,linewidth=1.0)
                plts.plot(the_lnlnX,the_lnlnY,c=color,linewidth=1.0,linestyle='dashed')
                

                #plts.hist(self.bears[i], 50, density=False, histtype='step', color=color)  #density = True if want integral of hist = 1
                #plts.hist(self.pears[i], 50, density=False, histtype='step', color=color,ls='--')  #probably not good to overlay with different bin sizes
            
                plts.set_xlabel(r'$ln(B/ \rho),lnB/ln \rho$')  
                plts.set_ylabel('PDF')
                
                # to get a percentage instead of count in the y axis for former .hist plotting usage
                y_vals = plts.get_yticks()
                plts.set_yticklabels(['{:.3f}'.format(x/len(the_lnX)) for x in y_vals])

            if 0: 
                # ALPHA FOR ALL PARTICLES OF ALL CORES PER FRAME
                if len(self.bears[i]) > 0: 
                    tmap2 = rainbow_map(len(self.bears[i]))
                    c2 = [tmap2(n) for n in range(len(self.bears[i]))]

                    plts.scatter(self.bears[i],self.pears[i],c=c2,alpha=0.2)   #WARNING: check what are bears and pears 
                    X = np.log10(self.bears[i])
                    Y = np.log10(self.pears[i])
                    X2 = np.linspace(X.min()+2,X.max()-3,num=len(X))  # EDIT, maybe take away the +2 and -3 
                    XX2 = 10 ** X2 
                   
                    # POLY-FITTING:
                    pfit = np.polyfit(X,Y,1)
                    other = np.poly1d(pfit)
                    alpha_fit = pfit[0]
                    print('alpha_fit ',alpha_fit)
                    B_o = pfit[1]  
                    YY = 10 ** (pfit[0]*X2 + pfit[1]) 
                    #plts.plot(XX2,YY,c='k',linewidth=1.0)
                    #xlims = 10e-3,10e7
                    #ylims = 10e-1,10e4
                    #BRho_tool.labelled(plts,xscale='log',yscale='log',xlabel=r'$\rho/\rho_{o}$',ylabel=r'$\mid B \mid (\mu G)$',\
                    #                    title=r'$\alpha_{fit} = %.3f$'%alpha_fit, xlim=xlims, ylim=ylims)

                    BRho_tool.labelled(plts,xscale='log',yscale=None,xlabel=r'$\rho$',ylabel=r'$log(B/Rho)$')#,\
                                        #title=r'$\alpha_{fit} = %.3f$'%alpha_fit, xlim=xlims, ylim=ylims)

            
            # EDIT RESP TO CALL LABELLED!!
            name_save = "TwoHistograms_%d_%d"%(i,num)  
            figs.savefig(name_save) 
            print('saved ',name_save)
            plt.close(figs)  
            figs, plts = plt.subplots(1,1) 
        
            # FOR WRITING STATS ON HISTOS
            #t = 'Cores: %d\n'%entries +  r'Mean $\beta = %.3f$'%betavg + '\n' + r'$\sigma = %.3f$'%betastd
            #ax1.text(0.56, 15, t, color='green', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u05 
            #ax1.text(-0.3, 3, t, color='green', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u101 
            #ax1.text(0.47, 13.5, t, color='blue', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u10 
            #ax1.text(0.43, 17.5, t, color='blue', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u102 
            #ax1.text(0.50, 12.5, t, color='m', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u11 
            #ax1.text(0.48,11, t, color='m', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u103  
            
    def boxes(self,num,name,tff_lab,ax): 
        Bears = {}
        index = []
        data = [] 
 
        for i,tff in enumerate(tff_lab): 
            #print('tff',tff)
            Bears[tff_lab[i]] = self.bears[i]
        for j, (key, val) in enumerate(Bears.items()):
            index.append(key)
            data.append(val)
        fig, ax1=plt.subplots(1,1)    

        meanpointprops = dict(marker='.', markeredgecolor='black',markerfacecolor='blue')  # was firebrick :), # marker options
        # showmeans = triangle green, meanline = solid orange line, by default, set meanline=True for dotted green of mean...
        #bparts = ax1.boxplot(data,showfliers=True,showmeans=True,meanline=True)  # bparts is now a dictionary; the originial box     
        bparts = ax1.boxplot(data,meanprops=meanpointprops,showmeans=True,meanline=False,showfliers=False)
        for i,val in enumerate(tff_lab):
            y = data[i] 
            x = np.random.normal(1+i, 0.04, size=len(y))
            ax1.plot(x, y, 'b.', alpha=0.2)
            if 1:
                y2 = data[i][-2]
                x2 = x[-2] 
                ax1.plot(x2,y2,'r.',alpha=1)  #alpha changes the opacity  
                #y3 = data[i][-1]
                #x3 = x[-1] 
                #ax1.plot(x3,y3,'k.',alpha=1)
                                          #CAREFUL because now I have changed the mean, see if it makes a huge difference :(,\
                                          #compare relative spacings with excel values
        #ax = df.boxplot()
        #ax.scatter(np.arange(df.shape[1])+1, df.loc[2000], color='r')
        #ax1.scatter(np.arange(ax1.shape[1])+1, bparts.loc[data[-1][-1]], color='r')  #EDIT


        # EXPERIMENTS
        if 0: #not tested yet
            q1 = pd.DataFrame(normal).quantile(0.25)[0]
            q3 = pd.DataFrame(normal).quantile(0.75)[0]
            iqr = q3 - q1 #Interquartile range
            fence_low = q1 - (1.5*iqr)
            fence_high = q3 + (1.5*iqr)
        if 0:
            means = [item.get_ydata()[0] for item in bparts['means']] 
            print(means,name) 
            if 0:
                alphaFile = open("alphaRecords.txt",'a')
                alphaFile.write("Sim %s meanAlphas %s \n"%(name,f'Means:{means}'))
                alphaFile.close()
        
        if 0:
            breakpoint()  #to try the dictionary
 
        #To compare boxplot with a zero, affects the position of the plots
        if num == 0: 
            ax1.plot([1,13],[0,0],c=[0.5]*4) 
            ax1.plot([1,13],[0.667,0.667],'m--',c=[0.5]*4)
            ax1.plot([1,13],[0.5,0.5],'--',c=[0.5]*4)
            ax1.plot([1,13],[0.4,0.4],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.07,0.15,0.23,0.30,0.38,0.45,0.53,0.61,0.68,0.75,0.82,0.89]#,0.95] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1))
            title = r'$\beta = 0.2$'
 
            therange = np.arange(0,1,0.075)
            #means = np.append(means,[np.nan]*1) 
            #ax.plot(therange,means,c='g',linestyle='dotted')
        if num == 1: 
            ax1.plot([1,12],[0,0],c=[0.5]*4) 
            ax1.plot([1,12],[0.667,0.667],'--',c=[0.5]*4)
            ax1.plot([1,12],[0.5,0.5],'--',c=[0.5]*4)
            ax1.plot([1,12],[0.4,0.4],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.08,0.16,0.25,0.33,0.41,0.50,0.58,0.66,0.74,0.82,0.90]#,0.97] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1))
            title = r'$\beta = 2.0$'

            therange = np.arange(0,1,0.075)
            #means = np.append(means,[np.nan]*2) 
            #ax.plot(therange,means,c='b',linestyle='dotted')
        if num == 2: 
            ax1.plot([1,11],[0,0],c=[0.5]*4) 
            ax1.plot([1,11],[0.667,0.667],'--',c=[0.5]*4)
            ax1.plot([1,11],[0.5,0.5],'--',c=[0.5]*4)
            ax1.plot([1,11],[0.4,0.4],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.08,0.17,0.26,0.35,0.44,0.52,0.60,0.69,0.77,0.86]#,0.91] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1)) 
            title = r'$\beta = 20$'

            therange = np.arange(0,1,0.075)
            #means = np.append(means,[np.nan]*3) 
            #ax.plot(therange,means,c='m',linestyle='dotted')

        ax1.set_xticklabels(tff_mod)       
        ylim = -0.5,1.0
        BRho_tool.labelled(ax1,xlabel=r'$t/t_{\rm{ff}}$',ylabel=r'$\alpha$',ylim=ylim,title=title) 
               
        fig.savefig('BearsBoxplot_%s'%name)
        print("BOX SAVED")
        plt.close(fig)
      

    def violins(self,num,name,tff_lab,ax):
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
         
        #vparts = ax1.violinplot(data,showmeans=True, showextrema=True, showmedians=True)
        #vparts['cmeans'].set_color('r')   
        #bparts = ax1.boxplot(data,showfliers=True,showmeans=True,meanline=True)

        meanpointprops = dict(marker='.', markeredgecolor='black',markerfacecolor='blue')  # was firebrick :), # marker options
        # showmeans = triangle green, meanline = solid orange line, by default, set meanline=True for dotted green of mean...
        #bparts = ax1.boxplot(data,showfliers=True,showmeans=True,meanline=True)  # bparts is now a dictionary; the originial box     
        bparts = ax1.boxplot(data,meanprops=meanpointprops,showmeans=True,meanline=False,showfliers=False)
        for i,val in enumerate(tff_lab):
            y = data[i] 
            x = np.random.normal(1+i, 0.04, size=len(y))
            ax1.plot(x, y, 'b.', alpha=0.2)
            y2 = data[i][-1]
            x2 = x[-1] 
            ax1.plot(x2,y2,'r.',alpha=1)

        # A TEST if bparts is used 
        if 0:
            means = [item.get_ydata()[0] for item in bparts['means']] 
            if 0:
                alphaFile = open("alphaRecords.txt",'a')
                alphaFile.write("Sim %s meanPears %s \n"%(name,f'Means:{means}'))
                alphaFile.close()
                
        #To compare violinplots with a zero, affects the position of the plots
        if num == 0:
            ax1.plot([1,13],[0,0],c=[0.5]*4)  
            ax1.plot([1,13],[1,1],'--',c=[0.5]*4)
            ax1.plot([1,13],[-1,-1],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.07,0.15,0.23,0.30,0.38,0.45,0.53,0.61,0.68,0.75,0.82,0.89]#,0.95] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1))
            title = r'$\beta = 0.2$'
            
            therange = np.arange(0,1,0.075)
            #means = np.append(means,[np.nan]*1) 
            #ax.plot(therange,means,c='g',linestyle='dotted')
        if num == 1:
            ax1.plot([1,12],[0,0],c=[0.5]*4) 
            ax1.plot([1,12],[1,1],'--',c=[0.5]*4)
            ax1.plot([1,12],[-1,-1],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.08,0.16,0.25,0.33,0.41,0.50,0.58,0.66,0.74,0.82,0.90]#,0.97] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1))
            title = r'$\beta = 2.0$'

            therange = np.arange(0,1,0.075)
            #means = np.append(means,[np.nan]*2) 
            #ax.plot(therange,means,c='b',linestyle='dotted')
        if num == 2:
            ax1.plot([1,11],[0,0],c=[0.5]*4) 
            ax1.plot([1,11],[1,1],'--',c=[0.5]*4)
            ax1.plot([1,11],[-1,-1],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.08,0.17,0.26,0.35,0.44,0.52,0.60,0.69,0.77,0.86]#,0.91] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1)) 
            title = r'$\beta = 20$'
       
            therange = np.arange(0,1,0.075)
            #means = np.append(means,[np.nan]*3) 
            #ax.plot(therange,means,c='m',linestyle='dotted')

        ax1.set_xticklabels(tff_mod)      
        ylim = -1.15,1.15 
        BRho_tool.labelled(ax1,xlabel=r'$t/t_{\rm{ff}}$',ylabel=r'$R$',ylim=ylim,title=title) 
         
        fig.savefig('PearsViolin_%s'%name)
        print("SAVED")
        plt.close(fig)



# - - - - - - - - - - - - - - - - - - - - 
# "MAIN"
# THREE SIMS AT ONCE

#import three_loopers_mountain_top as TLM
#import three_loopers_tenfour as TLTF
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


def axisforbox(theAx=None):
    # TO PLOT FIGURES OF ALL THREE SIMS AT ONCE, comment if overlaying with other py file
    #fig, ax1=plt.subplots(1,1)    
    for nt,tool in enumerate([BRho_tool1,BRho_tool2,BRho_tool3]):

        # TYPE OF PLOT: 'scatter_plot' OR 'frame_scatters' 
        # OR 'box_plot' OR 'vio_plot'? OR 'rms_plot' OR 'histogram'
        which_plot = 'vio_plot' 
        # ALL TIME: 'all_time', OR PER FRAME: 'per_frame'?
        which_time = 'per_frame'

        # GLOBAL TFF 
        G = 1620/(4*np.pi)
        rho_mean = 1
        t_ff = np.sqrt(3*np.pi/(32*G*rho_mean)) 

        # TFF PERCENTAGE & DESIRED FRAMES, DESIRED CORES
        # fix "core_list" respectively! 
        thtr = tool.this_looper.tr
        if nt == 0: 
            tff_p = thtr.times[:-1]/t_ff
            #tff_p = thtr.times[:]/t_ff
            frames = [1,3,6,8,10,12]   
            core_list = [27]
            corenum = 27
        if nt == 1:
            tff_p = thtr.times[:-1]/t_ff 
            #tff_p = thtr.times[:]/t_ff
            frames = [1,3,5,7,9,11] 
            core_list = [32]
            corenum = 32
        if nt == 2:
            tff_p = thtr.times[:-1]/t_ff 
            #tff_p = thtr.times[:]/t_ff
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
            print('here!')
            axnum,figu = tool.run(simnames[nt],thtr,which_time,which_plot,tff_labels,
                 fig=None,ax1=None,ax2= None,ax3= None,ax4= None,ax5= None,ax7= None,ax8= None,
                 lplots= None,core_list=None,core_val=corenum,simframes= None) 

        # check to see if these are now missing tool.run
        if 0:        
            # CALL HISTOGRAM FOR ALL CORES FOR ALL TIME or as in box-vios
            fig, ax = plt.subplots(1,1) 
            tool.histograms(nt,fig,ax)
        if 0: 
            # CALL BOXPLOTS  
            tool.boxes(nt,simnames[nt],tff_labels,theAx)
        if 1: 
            # CALL VIOLINPLOTS              
            tool.violins(nt,simnames[nt],tff_labels,theAx)
        if 0: 
            # TWO PANEL SCATTER PLOTS!
            tool.twopanels(nt, simnames[nt],axnum,figu,corenum) # EDIT

    # what was this for...
    # UNCOMMENT OUT WHEN USING THIS FILE TO RETURN A FIG
    #fig.savefig('BearsBoxplot_halphas')
    #print("BOXalphas SAVED")
    #plt.close(fig)


# COMMENT OUT WHEN USING THIS FILE AS AN IMPORT TO A PY FILE
axisforbox()
print("GOOD-BYE")



