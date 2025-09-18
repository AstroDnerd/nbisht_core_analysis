
import matplotlib
matplotlib.use('Agg')
import yt
import math
import matplotlib.pyplot as plt
import numpy as np
import pdb
nar = np.array

from importlib import reload

import looper
reload(looper)
import trackage
reload(trackage)
import tracks_read_write as trw
reload(trw)
from davetools import *
plt.close('all')
class mini_scrubber():
    def __init__(self,trk,core_id):
        self.trk=trk
        self.scrub(core_id)
        self.axis=0
                
    def scrub(self,core_id, axis=0):
        self.raw_x = self.trk.c([core_id],'x')
        self.raw_y = self.trk.c([core_id],'y')
        self.raw_z = self.trk.c([core_id],'z')
        #this_x=raw_x
        #this_y=raw_y
        if 1:
            #do the shift
            self.this_x = shift_down(self.raw_x)
            self.this_y = shift_down(self.raw_y)
            self.this_z = shift_down(self.raw_z)
        else:
            #don't actuall shift
            self.this_x = self.raw_x+0
            self.this_y = self.raw_y+0
            self.this_z = self.raw_z+0
        self.mean_x = np.mean(self.this_x,axis=0)
        self.mean_y = np.mean(self.this_y,axis=0)
        self.mean_z = np.mean(self.this_z,axis=0)
        self.nparticles,self.ntimes=self.this_x.shape
        self.meanx2 = np.tile(self.mean_x,(self.raw_x.shape[0],1))
        self.meany2 = np.tile(self.mean_y,(self.raw_x.shape[0],1))
        self.meanz2 = np.tile(self.mean_z,(self.raw_z.shape[0],1))
        self.r2 = (self.this_x-self.meanx2)**2+\
                  (self.this_y-self.meany2)**2+\
                  (self.this_z-self.meanz2)**2
        self.r=np.sqrt(self.r2)
        self.rmax = np.max(self.r,axis=0)
        self.max_track = np.where( self.r[:,0] == self.rmax[0])
        self.rmax_fat=np.tile(self.rmax,(self.raw_x.shape[0],1))
        self.rms = np.sqrt( np.mean(self.r2,axis=0))
        self.axis = axis
        if self.axis == 0:
            self.this_h = self.this_y
            self.h_label='y'
            self.this_v = self.this_z
            self.v_label='z'
        elif self.axis == 1:
            self.this_h = self.this_z
            self.h_label='z'
            self.this_v = self.this_x
            self.v_label='x'
        if self.axis == 2:
            self.this_h = self.this_x
            self.h_label='x'
            self.this_v = self.this_y
            self.v_label='y'
if 0:
    #newer read.
    if 'this_looper' not in dir():
        directory = '/home/dcollins/scratch/u05-r4-l4-128'
        this_looper=looper.core_looper(directory=directory)
        #file_list=glob.glob('/home/dcollins/scratch/Paper19/track_sixteen/*h5')
        file_list=glob.glob('/home/dcollins/scratch/Paper19/track_three/*h5')
        #file_list=glob.glob('/home/dcollins/scratch/Paper19/track_sixteen_good/*h5')

        for fname in file_list:

            trw.load_loop(this_looper,fname)
        thtr = this_looper.tr
        
import copy
def shift_up(pos,point):
    out = np.sort(copy.copy(pos),axis=1)
    out[ out < point ] = out[ out < point]+1
    return out
bork=0
delta=0
bo=0
smo=0
def shift_down(pos):
    #shift based on the time history: if a particle jumps more than half the box,
    #move it.
    global bork
    global delta
    global bo
    global smo
    sign = -1
    out = copy.copy(pos)
    bork=out[:,-1]+0
    #delta = sign*(out[:,1:]-out[:,:-1])
    #shape = delta.shape
    #bork=delta.max(axis=1)
    #bork = np.tile(bork,(shape[1],1)).transpose()
    bo = out+0#np.logical_and(delta <= bork, delta > 0.5)

    #out[:,:-1][bo] -= 1
   

    #distance_from_final = np.abs(out- np.tile(out[:,-1], (out.shape[1],1)).transpose())
    mean_pos = np.median(out[:,-1])
    print("mean_pos",mean_pos)
    distance_from_final =        out- mean_pos
    #ft = np.abs(distance_from_final[:,:-1]) > 0.5
    ft = np.abs(distance_from_final) > 0.5
    smo=ft
    out[ft] -=  1*np.sign(distance_from_final[ft])



    #out[ out > point ] = out[ out > point]-1
    return out#,delta
#moo = this_x[336:337,:]
#moo=this_x
#print(moo)
#cloo=shift_down(moo)
##plt.clf()
#plt.plot(moo[:,0],'k:')
#plt.plot(moo[:,1],'k-')
#plt.plot(moo[:,2],'k')
#plt.plot(cloo[:,-1],'g')
#plt.plot(cloo[:,0])
#plt.plot(cloo[:,1])
plt.savefig('image_tracks/dbg2.png')
if 0:
    #first read 
    if 'this_looper' not in dir():
        #directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
        directory = '/home/dcollins/scratch/u05-r4-l4-128'
        this_looper = looper.core_looper(directory= directory,
                                         sim_name = 'u05',
                                         out_prefix = 'test',
                                         target_frame = 125,
                                         frame_list = list(range(10,130,10)) + [125],
                                         core_list =  [ 120, 10, 308, 41, 44, 110, 1, 65], 
                                         fields_from_grid=['x','y','z']
                                      )
    if this_looper.tr is None:
        this_looper.tr = trackage.track_manager(this_looper)
        this_looper.tr.read('cores_many1.h5')
        thtr=this_looper.tr

#this_raw_x=raw_x[549,:]
#this_raw_y=raw_y[549,:]
#b=shift_down(this_raw_x)
##b=shift_down(this_y[0:1,:])
    ls=all_cores
#ls = [202]
#ls=[41]
#ls = [10]
#file_list=glob.glob('/home/dcollins/scratch/Paper19/track_three/*h5')
#file_list=glob.glob('/home/dcollins/scratch/Paper19/track_sixteen_good/*h5')
file_list=glob.glob('/home/dcollins/scratch/Paper19/track_sixteen_full/*h5')
file_list=glob.glob('/home/dcollins/scratch/Paper19/particle_error/particle_error_test_c0031_threeframes.h5')
file_list=glob.glob('/home/dcollins/scratch/Paper19/particle_error/track_indfix_sixteenframe_core_0031.h5')
directory = '/home/dcollins/scratch/u05-r4-l4-128'
core_31_baddies=nar([3192, 3207, 3283, 3327, 3390, 3444, 3458])
for nfile,fname in enumerate(file_list) :#[:3])
    #0164.h5
    t1 = fname.split("/")[-1]
    #l = len("track_three_to_test_core_")
    #l = len("track_sixteen_frames_core_")
    #this_cor = int(t1[l:l+4]) #[fname.index('_'):]
    this_cor=31
    if this_cor not in  [31]:
        continue
    print(this_cor)
    this_looper=looper.core_looper(directory=directory)
    trw.load_loop(this_looper,fname)
    thtr = this_looper.tr
    all_cores = np.unique(thtr.core_ids)
    ls=all_cores
    rm = rainbow_map(len(all_cores))
    plt.clf()
    plt.plot([0,1,1,0,0],[0,0,1,1,0])

    if 1:
        #time plots
        asort =  np.argsort(thtr.times)
        xxx = 1./128/np.logspace(0,4,5,base=2)
        fig_rt2,ax_rt2=plt.subplots(1,1)
        for nc,core_id in enumerate(ls):
            ms = mini_scrubber(thtr,core_id)
            density = thtr.c([core_id],'density')
            tmap=rainbow_map(ms.ntimes)
            for sub in range(20):
                plt.clf()
                for npart in range(sub*300,(sub+1)*300): #range(sub,ms.nparticles,20):
                    plt.plot( thtr.times[asort], density[npart,asort],c='k',linestyle=':',marker='*')
                outname = 'image_tracks/rho_t_sub2_%03d_c%04d.png'%(sub,core_id)
                plt.yscale('log')
                plt.savefig(outname)
                print('saved '+outname)

    if 0:
        #radial plots.
        asort =  np.argsort(thtr.times)
        xxx = 1./128/np.logspace(0,4,5,base=2)
        for nc,core_id in enumerate(ls):
            ms = mini_scrubber(thtr,core_id)
            density = thtr.c([core_id],'density')
            tmap=rainbow_map(ms.ntimes)
            if 1:
                #density plots
                tmap=rainbow_map(ms.ntimes)
                plt.clf()
                for xxxx in xxx:
                    plt.plot([xxxx,xxxx],[0.01,4e6],c=[0.5]*4)
                for it,nt in enumerate(asort):
                    plt.scatter(ms.r[:,nt], density[:,nt],c=[tmap(it)]*ms.r.shape[0],
                               label = "%0.4f"%thtr.times[nt])
                powerline(plt,1e-3,1e-1,5e4,-2)
                plt.xscale('symlog',linthreshx=xxx.min()/0.5);
                #plt.xlim(1e-5,0.5)
                plt.xlim(0,0.5)
                plt.yscale('log')
                plt.ylim(0.01,4e6)
                plt.legend(loc=1)
                outname = 'image_tracks/rho_t_c%04d.png'%core_id
                plt.savefig(outname)
                print('saved '+outname)

    
            #fig_many.savefig(outname)

    if 0:
        fig_rt, ax_rt = plt.subplots(1,1)
        ax0=ax_rt
        asort =  np.argsort(thtr.times)
        for nc,core_id in enumerate(ls):
            ms=mini_scrubber(thtr,core_id)
            tmap=rainbow_map(ms.ntimes)
            for num_part in range(5):#ms.nparticles):
                ax_rt.plot(thtr.times[asort],ms.r[num_part,asort])
            ax_rt.set_title('core %d'%core_id)
            outname = 'image_tracks/r_vs_time_c%04d.png'%core_id
            fig_rt.savefig(outname)
            print('saved '+outname)

    if 0:
        #Make images with circles.
        fig_rhist,ax_rhist=plt.subplots(1,1)
        fig,ax=plt.subplots(1,1,figsize=(8,8))
        fig_rmst,ax_rmst=plt.subplots(1,1,figsize=(8,8))
        ax.plot([0,1,1,0,0],[0,0,1,1,0])
        for nc,core_id in enumerate(ls):
            ax.clear()
            print('plot core %d'%core_id)
            n_for_color = int(np.where( all_cores == core_id)[0])
            ms = mini_scrubber(thtr,core_id)
            ax_rmst.plot(thtr.times, ms.rms,c=rm(n_for_color))
            ok = np.where( ms.r==ms.rmax_fat)
            circle_max = plt.Circle( (ms.mean_x[0],ms.mean_y[0]), ms.rmax[0], 
                                    color=rm(n_for_color),fill=False)
            ax.add_artist(circle_max)
            circle_rms = plt.Circle( (ms.mean_x[0],ms.mean_y[0]), ms.rms[0], 
                                    color=rm(n_for_color),fill=False)
            ax.add_artist(circle_rms)
            theseparts = np.arange(0,ms.nparticles,dtype='int')
            theseparts=core_31_baddies[0:1]
            rr=rainbow_map(len(theseparts))
            ax.plot([0,1,1,0,0],[0,0,1,1,0])
            for npp,npart in enumerate(theseparts):
                lab=None
                this_color = rm(n_for_color)
                if npart==0:
                    lab = 'c %d'%core_id
                    #print('wtf',lab)
                    this_color = rm(n_for_color)
                if 0:
                    this_color = rr(npp)
                ax.plot( ms.this_x[npart,:], ms.this_y[npart,:])#,c=this_color,label=lab)
            delta = 0.5
            ax.set_xlim(-delta,1+delta); ax.set_ylim(-delta,1+delta)
            outname = 'image_tracks/image_c%04d.png'%core_id
            fig.savefig(outname)
            print("save "+outname)

            tmap=rainbow_map(ms.ntimes)
            ax_rhist.clear()
            for nt in range(ms.ntimes):
                ax_rhist.hist(ms.r[:,nt],histtype='step',color=tmap(nt))
            ax_rhist.set_title('core %d'%core_id)
            outname = 'image_tracks/rhist_c%04d.png'%core_id
            fig_rhist.savefig(outname)
            print('saved '+outname)

            #lab = 'c %d'%core_id
            ##ax.plot(mean_x , mean_y,c=rm(nc),label=lab)
            #ax.plot(this_x-meanx2 ,this_y-meany2,c=rm(nc),label=lab)

        ax.set_xlim(-1,2); ax.set_ylim(-1,2)
        #ax.set_xlim(-delta,1+delta); ax.set_ylim(-delta,1+delta)
        ax.legend(loc=0)
        outname = 'image_tracks/core_rel.pdf'
        fig.savefig(outname)
        print('saved '+outname)
        outname = 'image_tracks/rrms.png'
        fig_rmst.savefig(outname)
        print('saved '+outname)
        plt.close('all')



