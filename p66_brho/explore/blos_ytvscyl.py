
'''
synthetic observations; version yt
print(ds.field_info['gas', 'magnetic_field_los'].get_source())

def _los_field(field, data):
    if data.has_field_parameter(f"bulk_{basename}"):
        fns = [(fc[0], f"relative_{fc[1]}") for fc in field_comps]
    else:
        fns = field_comps
    ax = data.get_field_parameter("axis")
    if is_sequence(ax):
    # Make sure this is a unit vector
        ax /= np.sqrt(np.dot(ax, ax))
        ret = data[fns[0]] * ax[0] + data[fns[1]] * ax[1] + data[fns[2]] * ax[2]
    elif ax in [0, 1, 2]:
        ret = data[fns[ax]]
    else:
        raise NeedsParameter(["axis"])
    return ret
'''

from cmath import log
from re import X
from turtle import width
from starter2 import *
import annotate_particles_4_cpy
reload(annotate_particles_4_cpy)
import multiplots 
reload(multiplots)
import davetools
reload(davetools)
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter
from matplotlib.ticker import PercentFormatter

np.set_printoptions(threshold=sys.maxsize)
from icecream import ic
# --- --- --- --- --- --- ---

class telescope(): 
    def __init__(self,the_loop):
        self.this_looper = the_loop
        self.cores_used = []

    def qtyRun(self,sim,core_list=None): 
        print('inside qtyRun')
        thtr = self.this_looper.tr

        # CORE_LIST
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores
        
        # FIELDS STORED
        self.synthRho = [np.zeros(len(core_list)) for x in range(3)]
        self.synthRho_mid = [np.zeros(len(core_list)) for x in range(3)]

        self.synthField = [np.zeros(len(core_list)) for x in range(3)]        
        self.synthField_mid = [np.zeros(len(core_list)) for x in range(3)]        

        self.synthField_yt_i = [np.zeros(len(core_list)) for x in range(3)]
        self.synthField_yt_ii = [np.zeros(len(core_list)) for x in range(3)]
        self.synthField_midyt_i = [np.zeros(len(core_list)) for x in range(3)]
        self.synthField_midyt_ii = [np.zeros(len(core_list)) for x in range(3)]

        # THE FINAL FRAME 
        the_frame = thtr.frames[-1:] 

        # CORE-LOOP
        for nc,core_id in enumerate(core_list):
            ds = self.this_looper.load(the_frame[0]) 
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            ms.particle_pos(core_id)
            self.ms = ms

            # PIECES FOR THE OBJECTS
            the_center = ms.mean_center[:,-1] 
            the_normal = [[1,0,0],[0,1,0],[0,0,1]]
            the_radius = 1/128
            the_area= np.pi * (the_radius**2) 

            # MAKE THE OBJECTS:
            xyz = [0,1,2]
            the_cyl = {}
            the_mid_cyl = {}

            for i in range(3):
                the_cyl[xyz[i]] = ds.disk(the_center,the_normal[i],the_radius,height=(1,'code_length'))
                the_mid_cyl[xyz[i]] = ds.disk(the_center,the_normal[i],the_radius,height=the_radius) 

                the_cyl[xyz[i]].set_field_parameter('axis',i)
                the_mid_cyl[xyz[i]].set_field_parameter('axis',i)

            # THE DENSITY & FIELD, ZONE METHOD & YT BLOS:  2D 
            B = ['magnetic_field_x','magnetic_field_y','magnetic_field_z']
            for j in range(3):
                self.synthRho[j][nc] = (the_cyl[j]['density'] * the_cyl[j]['cell_volume']).sum()/the_area
                self.synthRho_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j]['cell_volume']).sum()/the_area

                self.synthField[j][nc] = (the_cyl[j]['density'] * the_cyl[j][B[j]] * the_cyl[j]['cell_volume']).sum()/the_cyl[j]['gas','cell_mass'].sum()
                self.synthField_yt_i[j][nc] = (the_cyl[j]['magnetic_field_los'] * the_cyl[j]['cell_volume']).sum()
                self.synthField_yt_ii[j][nc] = (the_cyl[j]['density'] * the_cyl[j]['magnetic_field_los'] * the_cyl[j]['cell_volume']).sum()/the_cyl[j]['gas','cell_mass'].sum()

                self.synthField_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j][B[j]] * the_mid_cyl[j]['cell_volume']).sum()/the_mid_cyl[j]['gas','cell_mass'].sum()
                self.synthField_midyt_i[j][nc] = (the_mid_cyl[j]['magnetic_field_los'] * the_mid_cyl[j]['cell_volume']).sum()
                self.synthField_midyt_ii[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j]['magnetic_field_los'] * the_mid_cyl[j]['cell_volume']).sum()/the_mid_cyl[j]['gas','cell_mass'].sum()

 
        # CACHE TO DISK
        if 0:
            hfivename = 'p66_brho/h5files/blos_ytvscyls_%s.h5'%(sim)
            Fptr = h5py.File(hfivename,'w')
            Fptr['synthrho'] = self.synthRho
            Fptr['synthrhomid'] = self.synthRho_mid

            Fptr['synthblos'] = self.synthField
            Fptr['synthblosmid'] = self.synthField_mid
            
            Fptr['synthblos_yt_i'] = self.synthField_yt_i
            Fptr['synthblos_yt_ii'] = self.synthField_yt_ii
            Fptr['synthblosmid_yt_i'] = self.synthField_midyt_i
            Fptr['synthblosmid_yt_ii'] = self.synthField_midyt_ii
            Fptr.close


# MAIN
import three_loopers_six as TL6
if 'clobber' not in dir():
    clobber=True
if 'scope1' not in dir() or clobber:
    scope1=telescope(TL6.loops['u601'])
if 'scope2' not in dir() or clobber:
    scope2=telescope(TL6.loops['u602'])
if 'scope3' not in dir() or clobber:
    scope3=telescope(TL6.loops['u603'])

simnames = ['u601','u602', 'u603']

for nt,tool in enumerate([scope1,scope2,scope3]): 
    # WHICH CORES & RUN IF NOT CACHED
    if 0:
        all_cores = np.unique(tool.this_looper.tr.core_ids)
        #core_list = all_cores[:2]  #DEBUG
        core_list = all_cores
        tool.qtyRun(nt,core_list=core_list) 
        
    #ONCE CACHED, PLOTS 
    if 1:
        hfivename = 'p66_brho/h5files/blos_ytvscyls_%s.h5'%(nt)
        Fptr = h5py.File(hfivename,'r')
        synth_rho = Fptr['synthrho'][()]
        synth_rhomid = Fptr['synthrhomid'][()]

        synth_blos = Fptr['synthblos'][()]
        synth_blosmid = Fptr['synthblosmid'][()]

        synth_blosyt_one = Fptr['synthblos_yt_i'][()]
        synth_blosyt_two = Fptr['synthblos_yt_ii'][()]
        synth_blosmidyt_one = Fptr['synthblosmid_yt_i'][()]
        synth_blosmidyt_two = Fptr['synthblosmid_yt_ii'][()]

        dxyz = ['x','y','z']
        color=['b','g','orange']
        for i in range(3): 

            if 1:
                plt.scatter(synth_blosmid[i],synth_blosmidyt_one[i],c=color[i],alpha=0.7)
                plt.axline((0, 0), slope=1, c='k', linewidth=0.8)
                #plt.xscale('log')
                #plt.yscale('log')
                plt.xlabel(r'$B_{los:x,y,z}$ SHORT')
                plt.ylabel(r'$B_{los:x,y,z}$ SHORT, YTM1') 
                #plt.ylabel(r'$B_{los:x,y,z}$ LONG, YTrho') 
                #outname ='Bpol_vs_Btol_longlog_%d_%s'%(i,simnames[nt])
                outname ='Blos_ytm1_vs_us_short_%s'%(simnames[nt])

            # SAVE ONE AXIS AT A TIME
            if 0:
                plt.savefig(outname)
                print('plotted_%d'%i)
                plt.close('all')


        # SAVING ALL AXES AT ONCE
        if 1:
            plt.savefig(outname)
            print('plotted_%d'%i)
            plt.close('all')


        # AND CLOSE H5 FILE
        Fptr.close()

