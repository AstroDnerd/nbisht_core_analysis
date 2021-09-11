
# <B> vs <n>
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

                #self.mean_field_comps[core_id].append((bb * cell_volume).sum()/cell_volume.sum())
                self.mean_field_comps_py[core_id].append(bb.mean())  

                #self.mean_rho[core_id].append((density * cell_volume).sum()/(cell_volume.sum()))  
                self.mean_rho_py[core_id].append(density.mean())  

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
        outname='BnTracks_py_pl_%s'%simnames[nt]  #also do BnTest_py, one field at a time


        fig, ax1=plt.subplots(1,1) 
        the_x = np.empty([0],dtype=float)
        the_y = np.empty([0],dtype=float)

        # MAKE THE FIELDS INTO A 2D ARRAY WE CAN PLOT
        for ncore,core_id in enumerate(tool.cores_used):
            this_rho = tool.mean_rho_py[core_id] 
            this_field = tool.mean_field_comps_py[core_id] 
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
           
            ax1.scatter(this_rho[0], this_field[0],c='b',marker='*')  
            ax1.scatter(this_rho[-1], this_field[-1],c='r',marker='*')  
            ax1.plot(this_rho,this_field,c=[0.5]*4)

        # PLOT THE POWER LAW 
        pfit = np.polyfit(the_x,the_y,1)
        alpha = pfit[0]
        Bavg_o = pfit[1]

        X = np.linspace(the_x.min(),the_x.max(),num=len(the_x))  #short: -2, +3   
        XX = 10 ** X
        Y = 10 ** (pfit[0]*X + pfit[1])                
                                                       
        ax1.plot(XX,Y,c='k',linewidth=1.0)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel(r'$<\rho>$')
        ax1.set_ylabel(r'$<B>$')
        ax1.set_xlim(10e-2,10e7)
        ax1.set_ylim(10e-1,10e3)
        ax1.set_title(r'$\alpha = %.3f$'%alpha)
        fig.savefig(outname)
        print("saved")
'''
        # PLOT A FEW OF THE TRACKS
        nc = len(tool.cores_used)
        take_a_few = ((nc-1)*np.random.random(10)).astype('int')
        for ncore,core_id in enumerate(nar(tool.cores_used)[take_a_few]):

            # Adjust respectively 
            ax1.plot(this_rho,this_field,c=[0.5]*4) 

        fig.savefig(outname)
        print("saved")
'''          
