
from starter1 import *
import matplotlib
matplotlib.use('Agg')
import yt
import math
import matplotlib.pyplot as plt
import numpy as np
import pdb
import looper
import trackage
import time_kludge
from importlib import reload
reload(looper)
reload(trackage)
reload(time_kludge)

def get_savefile_type(savefile):
    #savefile can be "full" or "trackage"
    #"full" stores everyting and is very slow to read.
    #"trackage" stores only the track manager, and is very fast.
    fptr = h5py.File(savefile,'r')
    try:
        if "savefile_type" in fptr:
            savefile_type = fptr['savefile_type'][()]
        else:
            savefile_type = "full"
    except: #try A/except B/finally C: do A, if it fails do B, do C even if it fails. 
        raise
    finally:
        fptr.close()
    return savefile_type
    
def save_loop(self,fname):
    fptr = h5py.File(fname,'w')
    fptr['savefile_type'] = 'full'
    #just a bunch of values
    try:
        for value in [ "current_frame",
                       "data_template",
                       "sim_name",     
                       "frame_list",   
                       "core_list",    
                       "directory",    
                       "target_frame", 
                       "out_prefix",   
                       "plot_directory",
                      ]:
            result = self.__dict__[value]
            if result is None:
                continue
            fptr.create_dataset(value,data=result)
        fptr.create_dataset('fields_from_grid',data=np.array( [np.string_(s) for s in self.fields_from_grid]))
        ti_grp=fptr.create_group('target_indices')
        for core in self.target_indices:
            ti_grp.create_dataset(str(core),data=self.target_indices[core])

        fptr.create_group('snaps')
        for frame in self.snaps:
            frame_group_name = 'frame %d'%frame
            grp=fptr['snaps'].create_group(frame_group_name)
            for core in self.snaps[frame]:
                snp=grp.create_group("core %d"%core)
                this_snap = self.snaps[frame][core]
                for val in ["core_id",
                            "time",
                            "frame",
                            "mask",
                            "pos",
                            "ind",
                            "vel",
                            "R_centroid",
                            "R_vec",
                            "R_mag",
                            "N_vec",
                            "V_bulk",
                            "V_rel",
                            "V_radial"]:
                    result = this_snap.__dict__[val]
                    if result is not None:
                        snp.create_dataset(val,data=result)
                fv=snp.create_group('field_values')
                for val in this_snap.field_values:
                    fv.create_dataset(val, data=this_snap.field_values[val])

        tr_gr=fptr.create_group('track_manager')
        for val in ['particle_ids','core_ids','frames','times']:
            tr_gr.create_dataset(val,data= self.tr.__dict__[val])
        tr_di = tr_gr.create_group('track_dict')
        for k in self.tr.track_dict:
            tr_di[k] = self.tr.track_dict[k]
        #tr_gr[] = self.core_ids
        #tr_gr[] = self.frames
        #tr_gr[] = self.times
    except:
        raise
    finally:
        fptr.close()

def load_loop(self,fname):
    fptr = h5py.File(fname,'r')
    #just a bunch of values
    try:
        for value in [ "current_frame",
                       "data_template",
                       "sim_name",     
                       "frame_list",   
                       "core_list",    
                       "target_frame", 
                       "out_prefix",   
                       "fields_from_grid",
                       "output_directory"
                      ]:
            if value in fptr:
                the_value = fptr[value][()]
                if value in ['frame_list','core_list']:
                    for iii in the_value:
                        if iii not in self.__dict__[value]:
                            self.__dict__[value].append(iii)
                else:
                    self.__dict__[value] = the_value
        if self.directory is None:
            self.directory = fptr['directory'][()]
        #ti_grp=fptr.create_group('target_indices')
        #for core in self.target_indices:
        #    ti_grp.create_dataset(str(core),data=self.target_indices[core])
        for core in fptr['target_indices']:
            self.target_indices[int(core)]=fptr['target_indices'][core][()]

        for frame_name in fptr['snaps']:

            frame_number = int(frame_name.split()[1])
            self.load(frame_number,dummy=True)

            for core_name in fptr['snaps'][frame_name]:
                core_number = int(core_name.split()[1])
                print('LOADING', frame_number, core_number)
                core_grp = fptr['snaps'][frame_name][core_name]
                self.snaps[frame_number][core_number]=self.make_snapshot(frame_number,core_number,
                                                                        dummy_ds=True)
                this_snap=self.snaps[frame_number][core_number]
                if 'time' in core_grp:
                    this_snap.time = core_grp['time'][()]
                else:
                    this_snap.time = time_kludge.d[frame_number]
                #this_snap.field_values={}
                for val in ["core_id",
                            "frame",
                            "mask",
                            "N_vec"]:
                    if val in core_grp:
                        this_snap.__dict__[val]=core_grp[val][()]
                for val in ["R_centroid",
                            "R_vec",
                            "pos",
                            "ind",
                            "vel",
                            "R_mag", ]:
                    if val in core_grp:
                        #thisval = this_snap.ds.arr(core_grp[val].value,'code_length')
                        thisval = yt.units.yt_array.YTArray(core_grp[val][()],'cm')
                        this_snap.__dict__[val]=thisval
                for val in ["V_bulk",
                            "V_rel",
                            "V_radial"]:
                    if val in core_grp:
                        #thisval = this_snap.ds.arr(core_grp[val].value,'code_velocity')
                        thisval = yt.units.yt_array.YTArray(core_grp[val][()],'cm/s')
                        this_snap.__dict__[val]=thisval
                for val in core_grp['field_values']:
                    this_snap.field_values[val]=core_grp['field_values'][val][()]
                if self.tr is None:
                    self.tr = trackage.track_manager(self)
                self.tr.ingest(this_snap)

    except:
        raise
    finally:
        fptr.close()


def save_loop_trackage_only(self,fname):
    fptr = h5py.File(fname,'w')
    fptr['savefile_type'] = 'trackage'
    try:  #the try-except-finally is for clean failures with the hdf5 file

        #
        # write meta data
        #
        for value in [ "current_frame",
                       "data_template",
                       "sim_name",     
                       "frame_list",   
                       "core_list",    
                       "directory",    
                       "target_frame", 
                       "out_prefix",   
                       "plot_directory",
                      ]:
            result = self.__dict__[value]
            if result is None:
                continue
            fptr.create_dataset(value,data=result)
        # write fields
        fptr.create_dataset('fields_from_grid',data=np.array( [np.string_(s) for s in self.fields_from_grid]))
        # write target_indices
        ti_grp=fptr.create_group('target_indices')
        for core in self.target_indices:
            ti_grp.create_dataset(str(core),data=self.target_indices[core])

        tr_gr=fptr.create_group('track_manager')
        self.tr.write(fptr=tr_gr)
    except:
        raise
    finally:
        fptr.close()

def load_trackage_only(self,fname):
    fptr = h5py.File(fname,'r')
    #just a bunch of values
    try:
        for value in [ "current_frame",
                       "data_template",
                       "sim_name",     
                       "frame_list",   
                       "core_list",    
                       "target_frame", 
                       "out_prefix",   
                       "fields_from_grid"
                       "plot_directory",
                      ]:
            if value in fptr:
                the_value = fptr[value][()]
                if value in ['frame_list','core_list']:
                    if value not in self.__dict__:
                        self.__dict__[value]=[]
                    for iii in the_value:
                        if iii not in self.__dict__[value]:
                            self.__dict__[value].append(iii)
                else:
                    self.__dict__[value] = the_value
        if self.directory is None:
            self.directory = fptr['directory'][()]
        #ti_grp=fptr.create_group('target_indices')
        #for core in self.target_indices:
        #    ti_grp.create_dataset(str(core),data=self.target_indices[core])
        for core in fptr['target_indices']:
            self.target_indices[int(core)]=fptr['target_indices'][core][()]
        if 'track_manager' not in fptr:
            print("No track manager in this file.")
            raise
        self.tr = trackage.track_manager(self, h5ptr=fptr['track_manager'])

    except:
        raise
    finally:
        fptr.close()

def check(lpr):
    f0=lpr.frame_list[1]
    c0=lpr.core_list[1]
    for frame in lpr.frame_list[1:]:
        for core in lpr.core_list:
            if core not in lpr.snaps[frame]:
                print("missing d%d c%d"%(frame,core))
                continue
            for field in lpr.snaps[f0][c0].field_values:
                if field == 'V_radial':
                    continue
                fv0 = lpr.snaps[f0][core].field_values[field]
                #this_looper.snaps[0][77].field_values['density'].size
                fv1 = lpr.snaps[frame][core].field_values[field]
                if fv0.size != fv1.size:
                    print("n %d c %d f %s z %d %d"%(frame,core,field, fv0.size, fv1.size))



