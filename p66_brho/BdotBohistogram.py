
#
from starter2 import *
import data_locations as dl
from collections import defaultdict

import davetools
reload(davetools)

plt.close('all')


class magfield_density_tool2():  # was mass_tool()
    def __init__(self,this_looper):
        self.this_looper=this_looper
        
        self.mean_rho=defaultdict(list)
        self.mean_field_comps=defaultdict(list)  
        self.mean_angle = defaultdict(list)
        self.mean_BdotBo = defaultdict(list) 

        self.cores_used=[] 


    def run(self,simnames,core_list=None):
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
       
        # NOTES
        # got rid of the core loop to avoid loops...AND to get info straight from thtr
        # miniscrubber best for locations around centroid
        # do: mask & thtr

        fig, ax1=plt.subplots(1,1)
        colormap = rainbow_map(len(thtr.frames)) 
        #bins = np.linspace(-1500,1500,50)  #-1000,1000,50
        bins = np.linspace(-1,1,50)  #-1000,1000,50
        for nf,frame in enumerate(thtr.frames):  
           
            # Computes the unique particle mask.
            # Code note: we should get this tool to figure out dx by itself.
            # For all cores at once this time: 
            nx = int(1/dx)
            ix =np.floor(thtr.c(core_list,'x')/dx)[:,nf]
            iy =np.floor(thtr.c(core_list,'y')/dx)[:,nf]
            iz =np.floor(thtr.c(core_list,'z')/dx)[:,nf]
            index = ix + nx*(iy * nx*iz)
            ar = np.argsort(index)
            rs = np.argsort(ar)
            isorted=index[ar]
            mask = np.ones_like(ix,dtype='bool')
            mask[1:] = isorted[1:]-isorted[:-1] != 0
            mask2 = mask[rs]
          

            density = thtr.c(core_list,'density')[mask2,nf]
            cell_volume = thtr.c(core_list,'cell_volume')[mask2,nf]
 
            bx = thtr.c(core_list,'magnetic_field_x')[mask2,nf]
            by = thtr.c(core_list,'magnetic_field_y')[mask2,nf]
            bz = thtr.c(core_list,'magnetic_field_z')[mask2,nf]                
            bb = np.sqrt(bx*bx+by*by+bz*bz) 

            # <B(t=0)>
            bx_o = thtr.c(core_list,'magnetic_field_x')[mask2,0]
            by_o = thtr.c(core_list,'magnetic_field_y')[mask2,0]
            bz_o = thtr.c(core_list,'magnetic_field_z')[mask2,0]
            bb_o = np.sqrt(bx_o*bx_o + by_o*by_o + bz_o*bz_o)
            
            BtdotBo = bx*bx_o + by*by_o + bz*bz_o 
            #MAX = BtdotBo.max()
            #MIN = BtdotBo.min()
            #print("MAX",MAX)
            #print("MIN",MIN)
            costheta = BtdotBo/(bb*bb_o)

            outname='costheta__histograms_%d_%s'%(frame,simnames)
            print("curious")
            
            #ax1.hist(BtdotBo, bins=bins, density=False, histtype='step', color=colormap(nf))
            ax1.hist(costheta, bins=bins, density=False, histtype='step', color=colormap(nf))
        axbonk(ax1,xlabel='cos(theta)',ylabel='n',yscale='log')  #logging in histo bins
        fig.savefig(outname)
        print(outname)



#import three_loopers_mountain_top as TLM
import three_loopers_tenfour as TLTF
if 'clobber' not in dir():
    clobber=True
if 'mag_den1' not in dir() or clobber:
    mag_den1=magfield_density_tool2(TLTF.loops['u401'])
    simname1 = 'u401'   
if 'mag_den2' not in dir() or clobber:
    mag_den2=magfield_density_tool2(TLTF.loops['u402'])
    simname2 = 'u402' 
if 'mag_den3' not in dir() or clobber:
    mag_den3=magfield_density_tool2(TLTF.loops['u403'])
    simname3 = 'u403'

simnames = [simname1, simname2, simname3]

if 1:
    for nt,tool in enumerate([mag_den1,mag_den2,mag_den3]):
        # SET UP THE VARIABLE
        G=1620./(4*np.pi)
        tff_global = np.sqrt(3*np.pi/(32*G*1))

        tool.run(simnames[nt])

