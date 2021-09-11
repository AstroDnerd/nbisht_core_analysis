
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

        self.brms = np.empty([0],dtype=float)
        self.vrms = np.empty([0],dtype=float)
        self.vrmsq_bmean = np.empty([0],dtype=float)

        self.pears = [] 
        self.bears = []
   

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


    def run(self,name,thtr,timing,typeplot,
            fig=None,ax1=None,ax2=None,ax3=None,ax4=None,ax5=None,ax7=None,ax8=None,
            lplots=None,core_list=None,simframes=None):

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
                         lplots=None,xx2=None,yy=None):
           
            if typeplot == 'scatter_plot': 
                ax1.scatter(density[mask,n_time],magfield[mask,n_time],c=c,label=thtr.times[n_time],s=0.1)  #edit thtr.times          
                ax1.plot(xx2,yy,c='k',linewidth=1.0) #c=c[0] for colors
                
                BRho_tool.labelled(ax1,xscale='log',yscale='log',xlabel=r'$\rho/\rho_{o}$',ylabel=r'$\mid B \mid (\mu G)$',
                                   xlim=rho_extents.minmax, ylim=magfield_extents.minmax,title=r'$\alpha_{b\rho} = %.3f$'%beta)

                if n_time == asort[-1]:
                    outname = 'brhotff_all_c%04d_%s'%(core_id,name)  
                    plt.savefig(outname)  # CAREFUL WITH FIG VS PLT
                    print("saved "+outname)
                    plt.close(fig)  


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
                    c2 = [tmap2(n) for n in range(len(self.betarr))] 

                    if name == 'u301':
                        #tmap2 = rainbow_map(len(self.betarr[2:])) 
                        #c2 = [tmap2(n) for n in range(len(self.betarr[2:]))] 

                        the_range = np.arange(0,1,0.075) 
                        ax3.scatter(the_range,self.betarr,c=c2)   #[2:],
                        ax3.plot(the_range,self.betarr,c='k',linewidth=1.0)  
                        ax3.set_title(r"$\beta = 0.2, \alpha_{B}$")
                    if name == 'u302':
                        the_range = np.arange(0.075,1,0.075) 
                        ax3.scatter(the_range,self.betarr,c=c2)  
                        ax3.plot(the_range,self.betarr,c='k',linewidth=1.0)  
                        ax3.set_title(r"$\beta = 2.0, \alpha_{B}$")
                    if name == 'u303':
                        the_range = np.arange(0.075,0.9,0.075) 
                        ax3.scatter(the_range,self.betarr,c=c2)  
                        ax3.plot(the_range,self.betarr,c='k',linewidth=1.0)   
                        ax3.set_title(r"$\beta = 20, \alpha_{B}$")
                    ax3.plot([0,1],[0,0],c=[0.5]*4) 
                    ax3.plot([0,1],[0.4,0.4],'--',c=[0.5]*4)
                 
                    tff_mod = [0.0,0.0,0.3,0.6,0.9]  #this should match multiplelocator
                    ax3.set_xticklabels(tff_mod)   
                    ax3.xaxis.set_major_locator(plt.MultipleLocator(0.3))

                    ax3.yaxis.tick_right()
                    ax3.set_xlabel(r'$t_{\rm{ff}}$')  #\rm{  },to avoid italization 
                    ax3.set_yscale('linear')
                    ax3.set_ylim(-0.5,0.5)
          
                    self.betarr = np.empty([0],dtype=float)
   
                    outname = 'brhotff_c%04d_%s'%(core_id,name)  
                    plt.savefig(outname)  # CAREFUL WITH FIG VS PLT
                    print("saved "+outname)
                    plt.close(fig)  #CAREFUL, this doesn't allow to do another desired core :(  


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
                if name == 'u301':
                    self.bear12 = np.append(self.bear12,self.betarr[12])    
                    self.bear13 = np.append(self.bear13,self.betarr[13])
                    #self.bear14 = np.append(self.bear14,self.betarr[14])  
                self.bears = [self.bear0,self.bear1,self.bear2,self.bear3,self.bear4,self.bear5,self.bear6,\
                              self.bear7,self.bear8,self.bear9,self.bear10,self.bear11,self.bear12,self.bear13]#,self.bear14]
                self.betarr = np.empty([0],dtype=float) 
 
            if typeplot == 'vio_plot' and n_time == asort[-1]:    
                #print('NTIME_in',n_time) 
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
                self.pear11 = np.append(self.pear11,self.pearsonr[11])  #final u303 pear contains nans 
                if name == 'u301':
                    self.pear12 = np.append(self.pear12,self.pearsonr[12])  #final u302 pear, but want to discard final frame   
                    self.pear13 = np.append(self.pear13,self.pearsonr[13]) 
                    #self.pear14 = np.append(self.pear14,self.pearsonr[14]) 

                self.pears = [self.pear0,self.pear1,self.pear2,self.pear3,self.pear4,self.pear5,self.pear6,\
                             self.pear7,self.pear8,self.pear9,self.pear10,self.pear11,self.pear12,self.pear13]#,self.pear14]
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
            cv = thtr.c([core_id],'cell_volume')  
 
            if name != 'u301': 
                #print("name",name)
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


            # SETUP FIGURES FOR EACH CORE HERE: 
            #if which_plot == 'scatter_plot' or 'rms_plot':
                # FIX?
                #fig, ax1=plt.subplots(1,1)
                #print('needs editing')

            # - - - - - - - - - - - - - - - - - - - - 
            # TIMELOOP:
            for n_count,n_time in enumerate(asort):       
                mask = ms.compute_unique_mask(core_id, dx=1./2048,frame=n_time)  
                c=tmap(n_count,mask.sum())  #Used to be: ms.nparticles
                if timing == 'per_frame':   
                    X = np.log10(density[mask,n_time]).flatten()  # [mask,n_time], OR [:,n_time][mask]
                    Y = np.log10(magfield[mask,n_time]).flatten()
 
                if timing == 'all_time':      
                    time=thtr.times[n_time]
                    if time == 0:
                        continue
                    X = np.log10(density).flatten()  # unsure how to account for all time with 'mask'
                    Y = np.log10(magfield).flatten()

                if timing == 'all_time':
                    X2 = np.linspace(X.min()+2,X.max()-3,num=len(X))  # FOR ALL TIME, ONE PANEL PLOTS 
                if timing == 'per_frame':
                    X2 = np.linspace(X.min(),X.max(),num=len(X)) # FOR PER FRAME, MULTIPLE PANEL PLOTS 
                XX2 = 10 ** X2 
               
                # POLY-FITTING:
                pfit = np.polyfit(X,Y,1) 
                other = np.poly1d(pfit)
                beta = pfit[0]
                B_o = pfit[1]  
                YY = 10 ** (pfit[0]*X2 + pfit[1]) 

                if timing == 'all_time': 
                    if n_time == asort[-1]:
                        self.betarr = np.append(self.betarr,beta)
                else:
                    self.betarr = np.append(self.betarr,beta) 
                # THE NEGATIVE OUTLIERS
                if beta < 0: 
                    self.time_stamp = np.append(self.time_stamp,n_time)

                # CALL PEARSON R
                if typeplot == 'vio_plot':   
                    pearsonR(X,Y)
                    if n_time == asort[-1]:
                        scatterplots()                     
                        break  #TEMPORARY FIX to get violin plots :-/
                       
                # CALL SCATTER: ...unfortunately "or" behaves strangely
                if typeplot == 'frame_scatters':
                    scatterplots(fig,ax1,ax2,ax3,ax4,ax5,ax7,ax8,lplots,XX2,YY)  
                if typeplot == 'scatter_plot':  
                    scatterplots(fig,ax1,ax2,ax3,ax4,ax5,ax7,ax8,lplots,XX2,YY)  
                if typeplot == 'box_plot':
                    scatterplots()


                # RMS EXPERIMENTS
                if typeplot == 'rms_plot':  
                    if name != 'u301':  
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

                            if name == 'u302':
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
                                
                            if name == 'u303':
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


    def histograms(self,num,figs,plts):
    # HISTOGRAM FOR ALL CORES FOR ALL TIME

        #BETA ENTRIES, AVERAGE & STD 
        betavg = np.mean(self.betarr)
        betastd = np.std(self.betarr)
        entries = len(self.betarr) 
      
        if nt == 0:
            color = 'r'
        if nt == 1:
            color = 'b'
        if nt == 2:
            color = 'g'
        plts.hist(self.betarr, 50, density=False, histtype='step', color=color)

        if nt == 2:
            plts.set_xlabel(r'$\alpha_{B}$') 
            plts.set_ylabel('PDF')
            y_vals = plts.get_yticks()
            plts.set_yticklabels(['{:.3f}'.format(x/len(self.betarr)) for x in y_vals])

            t = 'Cores: %d\n'%entries +  r'Mean $\beta = %.3f$'%betavg + '\n' + r'$\sigma = %.3f$'%betastd
            #ax1.text(0.56, 15, t, color='green', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u05 
            #ax1.text(-0.3, 3, t, color='green', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u101 
            #ax1.text(0.47, 13.5, t, color='blue', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u10 
            #ax1.text(0.43, 17.5, t, color='blue', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u102 
            #ax1.text(0.50, 12.5, t, color='m', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u11 
            #ax1.text(0.48,11, t, color='m', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u103  
        
            # COULD CALL LABELLED!
            name_save = "BetaHistogramTff"  
            figs.savefig(name_save) 
            print('saved ',name_save)
            plt.close(figs)  


    def boxes(self,num,name,tff_lab): 
        Bears = {}
        index = []
        data = [] 
        for i,num in enumerate(tff_lab):              
            print(i,num)
            Bears[tff_lab[i]] = self.bears[i]
        for j, (key, val) in enumerate(Bears.items()):
            index.append(key)
            data.append(val)

        fig, ax1=plt.subplots(1,1)    
 
        bparts = ax1.boxplot(data,showfliers=True,showmeans=True,meanline=True)  #bparts is now a dictionary      

        # EXPERIMENTS
        if 0: #not tested yet
            q1 = pd.DataFrame(normal).quantile(0.25)[0]
            q3 = pd.DataFrame(normal).quantile(0.75)[0]
            iqr = q3 - q1 #Interquartile range
            fence_low = q1 - (1.5*iqr)
            fence_high = q3 + (1.5*iqr)
        if 0:
            print("Flier values of boxplot frame")
            flies = bparts['fliers'][0].get_ydata()
            print(flies)
        if 0:
            breakpoint()  #to try the dictionary

        #To compare boxplot with a zero, affects the position of the plots
        if nt == 0:
            ax1.plot([1,13],[0,0],c=[0.5]*4) 
            ax1.plot([1,13],[0.667,0.667],'m--',c=[0.5]*4)
            ax1.plot([1,13],[0.5,0.5],'--',c=[0.5]*4)
            ax1.plot([1,13],[0.4,0.4],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.07,0.15,0.23,0.30,0.38,0.45,0.53,0.61,0.68,0.75,0.82,0.89]#,0.95] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1))
            title = r'$\beta = 0.2$'
        if nt == 1:
            ax1.plot([1,12],[0,0],c=[0.5]*4) 
            ax1.plot([1,12],[0.667,0.667],'--',c=[0.5]*4)
            ax1.plot([1,12],[0.5,0.5],'--',c=[0.5]*4)
            ax1.plot([1,12],[0.4,0.4],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.08,0.16,0.25,0.33,0.41,0.50,0.58,0.66,0.74,0.82,0.90]#,0.97] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1))
            title = r'$\beta = 2.0$'
        if nt == 2:
            ax1.plot([1,11],[0,0],c=[0.5]*4) 
            ax1.plot([1,11],[0.667,0.667],'--',c=[0.5]*4)
            ax1.plot([1,11],[0.5,0.5],'--',c=[0.5]*4)
            ax1.plot([1,11],[0.4,0.4],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.08,0.17,0.26,0.35,0.44,0.52,0.60,0.69,0.77,0.86]#,0.91] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1)) 
            title = r'$\beta = 20$'

        ax1.set_xticklabels(tff_mod)      
        ylim = -1.0,1.0
        BRho_tool.labelled(ax1,xlabel=r'$t_{\rm{ff}}$',ylabel=r'$\alpha_B$',ylim=ylim,title=title) 
               
        fig.savefig('BearsBoxplot_%s'%name)
        print("BOX SAVED")
        plt.close(fig)
      

    def violins(self,num,name,tff_lab):
    # NOTE: for PearsonR, nan values are not accepted 
        Pears = {}
        index = []
        data = []  
        for i,num in enumerate(tff_lab):
            Pears[tff_lab[i]] = self.pears[i] 
        for j, (key, val) in enumerate(Pears.items()):
            index.append(key)
            data.append(val)

        fig, ax1=plt.subplots(1,1)    
 
        vparts = ax1.violinplot(data,showmeans=True, showextrema=True, showmedians=True)
        vparts['cmeans'].set_color('r')   

        #To compare violinplots with a zero, affects the position of the plots
        if nt == 0:
            ax1.plot([1,13],[0,0],c=[0.5]*4)  
            ax1.plot([1,13],[1,1],'--',c=[0.5]*4)
            ax1.plot([1,13],[-1,-1],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.07,0.15,0.23,0.30,0.38,0.45,0.53,0.61,0.68,0.75,0.82,0.89]#,0.95] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1))
            title = r'$\beta = 0.2$'
        if nt == 1:
            ax1.plot([1,12],[0,0],c=[0.5]*4) 
            ax1.plot([1,12],[1,1],'--',c=[0.5]*4)
            ax1.plot([1,12],[-1,-1],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.08,0.16,0.25,0.33,0.41,0.50,0.58,0.66,0.74,0.82,0.90]#,0.97] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1))
            title = r'$\beta = 2.0$'
        if nt == 2:
            ax1.plot([1,11],[0,0],c=[0.5]*4) 
            ax1.plot([1,11],[1,1],'--',c=[0.5]*4)
            ax1.plot([1,11],[-1,-1],'--',c=[0.5]*4)
            tff_mod = [0.0,0.0,0.08,0.17,0.26,0.35,0.44,0.52,0.60,0.69,0.77,0.86]#,0.91] 
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1)) 
            title = r'$\beta = 20$'
 
        ax1.set_xticklabels(tff_mod)      
        ylim = -1.15,1.15
        BRho_tool.labelled(ax1,xlabel=r'$t_{\rm{ff}}$',ylabel=r'$r_{B\rho}$',ylim=ylim,title=title) 
         
        fig.savefig('PearsViolinTen_%s'%name)
        print("SAVED")
        plt.close(fig)



# - - - - - - - - - - - - - - - - - - - - 
# "MAIN"
# THREE SIMS AT ONCE

import three_loopers_mountain_top as TLM
if 'clobber' not in dir():
    clobber=True

if 'BRho_tool1' not in dir() or clobber: 
    BRho_tool1=BRho_tool(TLM.loops['u301'])
    simname1 = 'u301' 

if 'BRho_tool2' not in dir() or clobber: 
    BRho_tool2=BRho_tool(TLM.loops['u302']) 
    simname2 = 'u302'
  
if 'BRho_tool3' not in dir() or clobber: 
    BRho_tool3=BRho_tool(TLM.loops['u303'])  
    simname3 = 'u303'

simnames = [simname1,simname2,simname3]
print("GREETINGS")


if 1:
    for nt,tool in enumerate([BRho_tool1,BRho_tool2,BRho_tool3]):

        # TYPE OF PLOT: 'scatter_plot' OR 'frame_scatters' 
        # OR 'box_plot' OR 'vio_plot'? OR 'rms_plot'
        which_plot = 'frame_scatters'
        # ALL TIME: 'all_time', OR PER FRAME: 'per_frame'?
        which_time = 'per_frame'

        # GLOBAL TFF 
        G = 1620/(4*np.pi)
        rho_mean = 1
        t_ff = np.sqrt(3*np.pi/(32*G*rho_mean)) 

        # TFF PERCENTAGE & DESIRED FRAMES, DESIRED CORES
        # FOR SCATTER PLOTS, comment 'core_list' if all cores: 
        thtr = tool.this_looper.tr
        if nt == 0: 
            tff_p = thtr.times[:-1]/t_ff
            frames = [1,3,5,8,10,12]   
            core_list = [275]  #this is not working, only does 1st core, following are blank.
        if nt == 1:
            tff_p = thtr.times[:-1]/t_ff 
            frames = [1,3,5,7,9,11] 
            core_list = [114]
        if nt == 2:
            tff_p = thtr.times[:-1]/t_ff 
            frames = [1,3,5,6,8,10] 
            core_list = [124]
        
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
                lplots = [0,ax1,0,ax2,0,ax4,0,0,ax5,0,ax7,0,ax8,0] 
            if nt == 1:
                lplots = [0,ax1,0,ax2,0,ax4,0,ax5,0,ax7,0,ax8,0] 
            if nt == 2:
                lplots = [0,ax1,0,ax2,0,ax4,ax5,0,ax7,0,ax8,0] 
     
            tool.run(simnames[nt],thtr,which_time,which_plot,
                 fig,ax1,ax2,ax3,ax4,ax5,ax7,ax8,
                 lplots,core_list,simframes) 
  
        if 0:  
            #if which_plot == 'scatter_plot' or 'rms_plot':  AVOID ORRRSSSSS
            #fig, ax1=plt.subplots(1,1)  #here or inside .run   

            # assign 'None' respectively 
            tool.run(simnames[nt],thtr,which_time,which_plot,
                 fig=None,ax1=None,ax2= None,ax3= None,ax4= None,ax5= None,ax7= None,ax8= None,
                 lplots= None,core_list=None,simframes= None) 

        # check to see if these are now missing tool.run
        if 0:        
            # CALL HISTOGRAM FOR ALL CORES FOR ALL TIME
            if nt == 0:
                fig, ax = plt.subplots(1,1) 
            tool.histograms(nt,fig,ax)
        if 0: 
            # CALL BOXPLOTS  
            tool.boxes(nt,simnames[nt],tff_labels)
        if 0: 
            # CALL VIOLINPLOTS              
            tool.violins(nt,simnames[nt],tff_labels)
    
    plt.close(fig)
print("GOOD-BYE")


# NEXT: how to close the figures for one core, but leave it open for the same in the same sim
'''
# TEMP CODE FOR INVESTIGATIONS
# in time loop, after X & Y in 'per_frame'
if n_time == 8:
    plt.hist(Y, 50, density=False, histtype='step', color='b')
    name_save = "MagfieldHistogramTff_colmask_%s"%name
    plt.savefig(name_save) 
    print('saved ',name_save)
    plt.close()  

    plt.hist(X, 50, density=False, histtype='step', color='g')
    name_save = "RhoHistogramTff_coldmask_%s"%name
    plt.savefig(name_save) 
    print('saved ',name_save)
    plt.close()  
'''
