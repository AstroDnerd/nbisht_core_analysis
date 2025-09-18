import matplotlib
matplotlib.use('Agg')
import yt
import math
import matplotlib.pyplot as plt
import numpy as np
nar = np.array
import time
import p19_tools as p19t
from p19_tools import *
reload(p19t)
def nothing(ds,sim,frame,ic, pos,vel,field_values, stop=False):
    return
def radial_velocity(ds,sim,frame,ic,leaf_indices,pos,vel,field_values, stop=False, other_args={}):
    m = field_values['density'].sum()
    shifted = p19t.shift_particles(ds,pos,shiftRight=False)


    R_centroid = ds.arr([0,0,0],'code_length')
    V_bulk = ds.arr([0,0,0],'code_velocity')
    for dim in range(3):
        R_centroid[dim] = (shifted[:,dim]*field_values['density']).sum()/m
        V_bulk[dim] = (vel[:,dim]*field_values['density']).sum()/m
    R_vec = shifted - R_centroid

    R_mag = (R_vec**2).sum(axis=1)**0.5
    N_vec = np.zeros_like(R_vec)
    for dim in range(3):
        N_vec[:,dim] = R_vec[:,dim]/R_mag
    V_relative = vel - V_bulk
    V_radial = (V_relative * N_vec).sum(axis=1)

    plt.scatter(R_mag, V_radial)
    for level in np.unique(ds.index.grid_levels):
        dx = ds.index.grids[0].dds[0]*0.5**level
        plt.plot([dx,dx],[V_radial.min(),V_radial.max()],c=[0.5]*4)
    plt.xlabel('R'); plt.ylabel(r'$V_{radial}$')
    if 'xscale' in other_args:
        plt.xscale(other_args['xscale'])
    if 'xlim' in other_args:
        xlim = other_args['xlim']
        plt.xlim(xlim)
        if R_mag.min() < xlim[0] or R_mag.max() > xlim[1]:
            print("X min warning: R_mag", R_mag.min(), R_mag.max())
    if 'ylim' in other_args:
        ylim = other_args['ylim']
        plt.ylim(ylim)
        if V_radial.min() < ylim[0] or V_radial.max() > ylim[1]:
            print("X min warning: V_radial", V_radial.min(), V_radial.max())
    outname = 'RadialVelocity_%s_c%04d_n%04d'%(sim,ic,frame)
    plt.savefig(outname); print(outname)
    if 'sphere_proj' in other_args:
        if other_args['sphere_proj']:
            for axis in [1]:
                proj = ds.proj('density',axis,center=R_centroid)
                pw = proj.to_pw(center = R_centroid, origin='domain')
                this_width = 2*R_mag.max().v
                this_width = 0.25
                pw.set_width((this_width,'code_length'))
                pw.set_cmap('density','gray')
                radius_from_clump = []
                clump_label = "c%04d_"%ic
                for nc,clump_number in enumerate(leaf_indices):
                #    clump_label += "%d_"%clump_number
                    pw.annotate_select_particles(1.0, col='r', indices=leaf_indices[clump_number])
                outname = '%s_no_particles_fc_%sn%04d'%(sim,clump_label,frame)
                #pw.annotate_grids()
                pw.annotate_sphere(R_centroid,0.05) #R_mag.max())
                print("SPHERE", R_centroid, R_mag.max())
                print pw.save(outname)




    if stop:
        pdb.set_trace()
