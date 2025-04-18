from starter2 import *
from yt.data_objects.level_sets import *
from collections import defaultdict
import weakref
import loop_tools
reload(loop_tools)
import tracks_read_write
reload(tracks_read_write)
verbose=True

def density_with_sinks(field, data):
    return data["gas", "density"] + (data["deposit", "all_density"])

def load_looper(fname, directory=None, mode_fname=None, density_placeholder='density'):
    import looper
    h5ptr=h5py.File(fname,'r')
    looper_version=1
    if 'looper_version' in h5ptr:
        looper_version=h5ptr['looper_version'][()]

    if looper_version==2:
        looper_main = core_looper2
    else:
        looper_main = looper.core_looper

    if directory is None:
        #It's actually not great to store this if you move machines.
        directory = h5ptr['directory'].asstr()[()]
    new_looper=looper_main(directory=directory, savefile_only_trackage=fname, density_placeholder = density_placeholder)
    new_looper.core_ids = np.sort(np.unique(new_looper.tr.core_ids))

    new_looper.mode_fname = mode_fname
    new_looper.read_mode()

    return new_looper



def count_particles(fname='tools_data/n_particles.txt'):
    fptr = open(fname,'r')
    lines=fptr.readlines()
    fptr.close()
    parts = np.zeros([len(lines),2])
    pd = {}
    for n,line in enumerate(lines):
        parts[n] = np.array(line.split(),dtype='int')
    pd = dict( zip( parts[:,0], parts[:,1]))
    return pd



def get_all_nonzero(fname='tools_data/n_particles.txt'):
    fptr = open(fname,'r')
    lines=fptr.readlines()
    fptr.close()
    parts = np.zeros([len(lines),2])
    for n,line in enumerate(lines):
        parts[n] = np.array(line.split(),dtype='int')
    all_nonzero = parts[:,0][ parts[:,1] >0]
    core_list = all_nonzero.astype('int')[::-1]
    return core_list
            
class core_looper2():
    """We want to analyze the core pre-images for a bunch of frames.
    *core_looper* keeps track of stuff (where the data is, how to loop, all the stuff)
    *snapshot* is for one core in one frame, and packs up the particle values.

    Usage:
        Make a core_looper with a 
         data_directory (where the stuff is)
         data_template  (how to get data from frames)
         out_prefix     (what to call your stuff)
         frame_list     (what frames to analyze)
         target_frame   (what frame the end states are in)
         
        and some extra stuff,
         fields_from_grid (in addition to density and cell volume)
         individual_particle_tracks (doesn't work),

        for example
         >>> u05_loop = looper.core_looper(directory='my_scratch',out_prefix='u05', frame_list=[0,10])

        You need to get the target indices with get_target_indices, this probably needs
        to be worked out better; right now pulling the peaks from an hdf5 file is best.
         >>> u05_loop.get_target_indices(h5name = 'u05_0125_peaklist.h5', core_list=[10])

        Then make some analysis loops, like the ones in run_loop.py
         
    """
    def __init__(self,savefile=None,
                 savefile_only_trackage=None,
                 directory="./", data_template = "%s/DD%04d/data%04d", sim_name='sim',
                 plot_directory="./plots_to_sort",
                 out_prefix="",
                 frame_list=None, core_list=None, target_frame=0,
                 fields_from_grid = None, 
                 individual_particle_tracks=False,
                 derived=None, do_shift=True, 
                 bad_particles=None,
                 density_placeholder = 'density'):
        #set defaults and/or arguments.
        self.looper_version = 2
        self.current_frame = None
        self.data_template = data_template
        self.sim_name      = sim_name
        if frame_list is None:
            frame_list=[]
        self.frame_list     = frame_list
        if core_list is None:
            core_list=[]
        self.core_list     = core_list
        self.directory     = directory
        self.target_frame  = target_frame
        self.out_prefix    = out_prefix
        self.plot_directory = plot_directory
        if fields_from_grid is None:
            fields_from_grid = []
        if density_placeholder == 'total_density':
            self.fields_from_grid = [YT_total_density, YT_cell_volume] + fields_from_grid
        else:
            self.fields_from_grid = [YT_density, YT_cell_volume] + fields_from_grid
        
        self.density_placeholder = density_placeholder

        #the track manager.
        self.tr = None

        #defaults for things to be set later
        self.target_indices = nar([])
        self.core_ids       = nar([])
        self.targets = None
        self.ds = None
        self.field_values=None
        self.snaps = {}
        #    defaultdict(whatev) is a dict, but makes a new (whatev) by default
        self.ds_list={}
        self.all_data={}
        if derived is None:
            derived=[]
        self.derived=derived

        self.shift = do_shift

        if bad_particles is None:
            bad_particles = defaultdict(list)
        self.bad_particles = bad_particles

        if savefile_only_trackage is not None:
            tracks_read_write.load_trackage_only(self,savefile_only_trackage )
        
        self.sink_core_dic = {}

        #read from save file.
        #if savefile is not None:
        #    fptr = open(savefile,'r')
        #    lines = fptr.readlines()
        #    fptr.close()
        #    for line in lines:
        #        line = line.strip() #clean off white space
        #        #treat each line as a command to run on 'self'
        #        exec("self.%s"%sline)

    def fix_frame_list(self):
        #Sometimes runs will have consecutive datasets with same time coz of restart
        frame_to_remove = []
        last_time = -1
        for frames in self.frame_list:
            with open(str(self.data_template)%(self.directory,frames,frames)) as datafile:
                for line in datafile:
                    if 'InitialTime' in line:
                        current_time = float(line.split()[-1])
                        if  current_time== last_time:
                            frame_to_remove.append(frames)
                        last_time = current_time
                        break
        print("Removing Frames!: ",frame_to_remove)
        self.frame_list = np.array(self.frame_list)[~np.isin(self.frame_list, frame_to_remove)]
        print("New Frame List: ",self.frame_list)

    def read_mode(self):
        import read_mode
        reload(read_mode)
        if self.mode_fname is None or not os.path.exists(self.mode_fname):
            print("No extant mode.  Please create.  Ignoring.")
            return
        if self.mode_fname is not None:
            package = read_mode.read( self.mode_fname)
            self.modes = package['modes']
            self.mode_dict=package['mode_dict']
            self.mode_label_dict=package['mode_label_dict']
            self.core_by_mode = package['core_by_mode']
            self.unique_modes = package['unique_modes']
            if np.abs( nar(sorted(package['core_ids'])) - nar(sorted(self.core_ids))).sum() > 0:
                print("Core ID problem.")
                raise
    def get_current_frame(self):
        if self.current_frame is None:
            self.current_frame = self.frame_list[0]
        return self.current_frame

    def load(self,frame=None,dummy=False,derived=None):
        """runs yt.load, and saves the ds so we don't have multiples for many cores."""
        if dummy:
            self.ds = None
            self.ds_list[frame]=None
            return None
        if frame is None:
            frame = self.get_current_frame()
        self.filename = str(self.data_template)%(self.directory,frame,frame)
        new_ds = True
        if frame in self.ds_list:
            if self.ds_list[frame] is not None:
                self.ds = self.ds_list[frame]
                new_ds = False
        if new_ds:
            self.ds = yt.load(self.filename)
            self.ds.add_field(("gas",self.density_placeholder),units="g/cm**3",function=density_with_sinks,sampling_type="cell")
            if derived is None:
                derived = self.derived
            for add_derived_field in derived:
                add_derived_field(self.ds)
            self.ds_list[frame] = self.ds
        if False:
            self.ds.force_periodicity(True)
        return self.ds
    def get_all_particles(self,frame):
        region = self.get_region(frame)
        self.all_particle_index =    region[YT_particle_index].astype('int64')
        self.all_particle_position = region[YT_particle_position]
        self.all_particle_velocity = region[YT_particle_velocity]

    def get_target_indices(self,target_frame=None,core_list=None,h5_name=None, peak_radius=1.5,
                          bad_particle_list=None):
        if target_frame is None:
            target_frame = self.target_frame
        if core_list is None and self.core_list is not None:
            core_list = self.core_list
        target_ds = self.load(target_frame)
        print("Stuff: ds %s target frame %s cores %s"%(str(target_ds), str(target_frame), str(core_list)))
        new_indices = loop_tools.get_leaf_indices(target_ds,h5_name = h5_name, 
                                     subset = core_list, peak_radius=peak_radius,
                                                 bad_particle_list=bad_particle_list)
        #for core_id in new_indices:
        #    self.target_indices[core_id] = new_indices[core_id]

        #make a single list of the indices,
        #and a list of core_ids for each particle.
        self.target_indices = np.concatenate([new_indices[core_id] for core_id in new_indices])
        self.core_ids       = np.concatenate([new_indices[core_id]*0+core_id for core_id in new_indices])

        #Sort the particles by index.
        args = np.argsort(self.target_indices)
        self.target_indices = self.target_indices[args]
        self.core_ids = self.core_ids[args]
        print('SORTING CORE IDS')
        

    def get_region(self,frame=None):
        if frame is not None or self.ds is None:
            ds = self.load(frame)
        region = self.ds.all_data()
        return region

    def make_snapshot(self,frame,dummy_ds=False):
        if frame in self.snaps:
            this_snap = self.snaps[frame]
        else:
            this_snap = snapshot(self,frame,dummy_ds=dummy_ds, density_placeholder = self.density_placeholder)
            self.snaps[frame] = this_snap # not a weak ref, needs to persist.weakref.proxy(this_snap)
        return this_snap

    def read_targets_only(self,target_fname):
        import mountain_top
        if self.targets is None:
            self.targets = {}
        h5ptr = h5py.File(target_fname,'r')
        for group_name in h5ptr:
            core_id = h5ptr[group_name]['peak_id'][()]
            #read from disk
            self.targets[core_id] = mountain_top.target_info(h5ptr=h5ptr[group_name])
        
    #reload(mountain_top)
    def read_targets(self,target_fname, bad_particle_fname_read):
        import mountain_top
        all_particles=np.array([])
        core_ids_by_particle=np.array([])
        if self.targets is None:
            self.targets = {}
        h5ptr = h5py.File(target_fname,'r')
        for group_name in h5ptr:
            core_id = h5ptr[group_name]['peak_id'][()]
            #read from disk
            self.targets[core_id] = mountain_top.target_info(h5ptr=h5ptr[group_name])
            these_targets= self.targets[core_id].particle_index
            self.target_indices = np.concatenate([self.target_indices, these_targets]).astype('int')
            this_id_array = np.array([core_id]*these_targets.size)
            self.core_ids = np.concatenate([self.core_ids, this_id_array]).astype('int')
            #self.target_indices = self.targets[core_id].particle_index

            #check for uniqueness
            if np.unique(self.target_indices).size - self.target_indices.size:
                arr, count_arr = np.unique(self.target_indices, return_counts=True)
                dup = arr[count_arr>1]
                bad_particles=bad_particle_fname_read
                bad_p = h5py.File(bad_particles, 'a')
                bad_p[str(core_id)] = dup
                bad_p.close()
                print("FATAL ERROR: %d repeated particles in core %d"%(len(dup),core_id))
                #pdb.set_trace()

        h5ptr.close()
        self.core_list = np.sort( np.unique( self.core_ids) )

        args = np.argsort(self.target_indices)
        self.target_indices = self.target_indices[args]
        self.core_ids = self.core_ids[args]
        print('SORTING CORE IDS')
        

    def verify_all_particles(self,frame):

        if self.bad_particles is None:
            self.bad_particles = defaultdict(list)

        all_target_particles = self.target_indices.astype('int64')
        if type(all_target_particles) == yt.units.yt_array:
            #strip off the units since we'll pass this to cython.
            all_target_particles = all_target_particles.v

        #check for repeated particles.  
        if np.abs(np.unique(all_target_particles).size - all_target_particles.size):
            print("Fatal error: repeated/missing particle.")
            raise

        #we'll need to know which cores have missing particles,
        #so make a list now.
        core_id_by_particle = self.core_ids

        #a mask of the particles we're looking for.  
        #records which ones are found.  Also makes things faster.
        mask_to_get=np.zeros(all_target_particles.shape,dtype='int32')

        #all indices in the simulation.  Sometimes misses particles.
        self.get_all_particles(frame)
        all_indices = self.all_particle_index.astype('int64')
        all_order = np.argsort(all_indices)
        all_return = np.argsort(all_order)
        sorted_all = all_indices[all_order]

        target_order = np.argsort( all_target_particles )
        sorted_targets = all_target_particles[target_order]
        return_targets = np.argsort(target_order)
        found_any, all_mask = particle_ops.mask_particles_sorted_t7(sorted_targets,sorted_all,mask_to_get)
        #mask = all_mask[all_return].astype('bool') #we don't actually need this here

        missing_particles = all_target_particles[ mask_to_get[return_targets] == 0]
        missing_cores     = core_id_by_particle[ mask_to_get[return_targets] == 0]

        for particle_id, core_id in zip(missing_particles,missing_cores):
            self.bad_particles[core_id].append( particle_id)

        bad_index = loop_tools.check_particles(self.ds_list[frame])
        for particle_id in bad_index:
            loc = np.where( all_target_particles == particle_id)[0]
            if len(loc) == 1:
                core_id = core_id_by_particle[loc[0]]
                self.bad_particles[core_id].append(particle_id)



    def remove_bad_particles(self):

        keepers = np.ones( self.target_indices.size, dtype='bool')
        for core_id in self.bad_particles:
            these_bad_particles = self.bad_particles[core_id]
            if core_id not in self.core_ids:
                continue
            for particle in self.bad_particles[core_id]:
                found_it =  np.where( self.target_indices == particle)
                keepers[found_it] = False
        self.target_indices = self.target_indices[keepers]
        self.core_ids = self.core_ids[keepers]

    def save_bad_particles(self,fname):
        h5ptr = h5py.File(fname,'w')
        for core_id in self.bad_particles:
            h5ptr[str(core_id)] = nar(self.bad_particles[core_id])
        h5ptr.close()
    def read_bad_particles(self,fname,core_hijack=None):
        h5ptr = h5py.File(fname,'r')
        if self.bad_particles is None:
            self.bad_particles = defaultdict(list)
        core_ids = []
        for group in h5ptr:
            core_ids.append(group)
        core_ids = np.sort(core_ids).astype('int')
        for core_id in core_ids:
            this_core_id = int(core_id)
            if core_hijack is not None:
                #this is a dirty kludge to get all particles to work.
                #
                this_core_id = core_hijack
            self.bad_particles[this_core_id] = h5ptr[str(core_id)][()]
        h5ptr.close()


            

    def get_tracks(self):
        if self.tr is None:
            self.tr = trackage.track_manager(self)
        for frame in self.frame_list:
            self.get_all_particles(frame)
            this_snapshot = self.make_snapshot(frame)
            if this_snapshot.R_centroid is None:
                this_snapshot.get_all_properties()
                #get_current_mask
                #get_pos_vel
                #get_particle_values_from_grid

            #this_snapshot.get_particle_values_from_grid()
            error=np.abs(np.sort(this_snapshot.ind) - this_snapshot.ind)
            if error.max() != 0:
                print("FATAL ERROR: particle order error")
                pdb.set_trace()
            self.tr.ingest(this_snapshot)

"""
        #this is not used.
        self.individual_particle_tracks=individual_particle_tracks

        #the track manager.
        self.tr = None

        #defaults for things to be set later
        self.target_indices = {}
        self.ds = None
        self.field_values=None
        self.snaps = defaultdict(dict) 
        #    defaultdict(whatev) is a dict, but makes a new (whatev) by default
        self.ds_list={}
"""



class snapshot():
    """For one core and one time, collect the particle positions and whatnot.
    """
    def __init__(self,loop,frame,dummy_ds=False, density_placeholder = 'density'):
        self.loop           = weakref.proxy(loop) #cyclic references are bad, weakref helps.
        self.target_indices = weakref.proxy(loop.target_indices)
        self.core_ids = weakref.proxy(loop.core_ids)
        self.dummy_ds=dummy_ds
        if dummy_ds:
            self.ds=None
            self.time = -1
        else:
            self.ds             = weakref.proxy(loop.load(frame,dummy=dummy_ds) )
            self.time = self.ds['InitialTime']
        self.frame          = frame

        #stubs for the quantities we'll compute.
        self.mask=None
        self.pos=None
        self.R_centroid=None
        self.field_values={}
        self.R_centroid  =None #(centroid weighted by grid quantities)
        self.R_vec       =None #(particle position relative to centroid)
        self.R_mag       =None #(magnitude of position)
        self.N_vec       =None #(normal vector)
        self.V_bulk      =None #(mean motion)
        self.V_rel       =None #(relative motion)
        self.V_radial    =None #(radial coordinate of velocity)
        
        self.dummy_ds = dummy_ds

        self.density_placeholder = density_placeholder

    def get_ds(self,frame=None):
        if frame is None:
            frame = self.frame
        self.ds             = weakref.proxy(self.loop.load(frame,dummy=False) )
        return self.ds

    def get_all_properties(self):
        """Run all the relevant analysis pieces."""
        print("get data frame %d"%(self.frame))
        self.get_current_mask()
        self.get_current_pos_vel()
        self.get_particle_values_from_grid()
        #self.compute_relative_coords()

    def get_region(self,region_stuff=None):
        """I will probably extend this to make more complex regions."""
        region = self.ds.all_data()
        return region

    def get_current_mask(self):
        """get the particle mask that relates particles for this core_id and this frame
        to the particles in the target_indices from the target_frame"""
        target_ids=self.target_indices.astype('int64')
        if type(target_ids) == yt.units.yt_array:
            target_ids = target_ids.v

        mask_to_get=np.zeros(target_ids.shape,dtype='int32')
        all_indices = self.loop.all_particle_index
        
        target_order = np.argsort( target_ids)
        sorted_targets = target_ids[target_order]
        all_order = np.argsort(all_indices)
        all_return = np.argsort(all_order)
        sorted_all = all_indices[all_order]
        found_any, all_mask = particle_ops.mask_particles_sorted_t7(sorted_targets,sorted_all,mask_to_get)

        if mask_to_get.sum() != mask_to_get.size:
            pdb.set_trace()
        self.mask = all_mask[all_return].astype('bool')
        return found_any, self.mask
        

    def get_current_pos_vel(self):

        if self.mask is None:
            found_any, mask = self.get_current_mask()
            if not found_any:
                raise
        else:
            mask = self.mask.astype('bool')
        region = self.get_region(self.frame)
        self.ind = self.loop.all_particle_index[mask]
        args=np.argsort(self.ind)
        self.ind=self.ind[args]
        self.pos = self.loop.all_particle_position[mask][args]
        self.vel = self.loop.all_particle_velocity[mask][args]


    def compute_relative_coords(self, density_placeholder = 'density'):
        """Compute and store 
        R_centroid (centroid weighted by grid quantities)
        R_vec      (particle position relative to centroid)
        R_mag      (magnitude of position)
        N_vec      (normal vector)
        V_bulk     (mean motion)
        V_rel      (relative motion)
        V_radial   (radial coordinate of velocity)
        Assumes that the following have been set properly:
        self.field_values
        self.pos
        self.vel
        self.ds
        """
        if len(self.field_values.keys())==0:
            self.get_particle_values_from_grid()

        m = self.field_values[density_placeholder].sum()
        if self.loop.shift:
            shifted = loop_tools.shift_particles(self.ds,self.pos,shiftRight=False)
        else:
            shifted = copy.copy(self.pos)
        if self.ds is not None:
            self.R_centroid = self.ds.arr([0,0,0],'code_length')
            self.V_bulk = self.ds.arr([0,0,0],'code_velocity')
        else:
            self.R_centroid = yt.units.yt_array([0,0,0],'cm')
            self.V_bulk =    yt.units.yt_array([0,0,0],'cm/s')
        centroid_tmp =  np.array([(shifted[:,dim]*self.field_values[density_placeholder]).sum()/m for dim in range(3)])
        vbulk_tmp = [ (self.vel[:,dim]*self.field_values[density_placeholder]).sum()/m for dim in range(3)]
        if self.ds is not None:
            self.R_centroid = self.ds.arr(centroid_tmp,'code_length')
            self.V_bulk = self.ds.arr(vbulk_tmp,'code_velocity')
        else:
            self.R_centroid = yt.units.yt_array(centroid_tmp,'cm')
            self.V_bulk =    yt.units.yt_array(vbulk_tmp,'cm/s')

        self.R_vec = shifted - self.R_centroid
        self.R_mag = (self.R_vec**2).sum(axis=1)**0.5
        self.N_vec = np.zeros_like(self.R_vec)
        for dim in range(3):
            self.N_vec[:,dim] = self.R_vec[:,dim]/self.R_mag
        self.V_relative = self.vel - self.V_bulk
        self.V_radial = (self.V_relative * self.N_vec).sum(axis=1)
        #self.field_values['V_radial']=self.V_radial
        #self.field_values['V_mag']=((self.V_relative*self.V_relative).sum())**0.5
        #self.field_values['R_mag']=self.R_mag
        #self.field_values['V_radial']=self.V_radial
        #self.field_values['V_relative']=self.V_relative


    def get_particle_values_from_grid(self, field_list=[]):
        """Get the simple nearest-sampled grid point from the particles in self.pos.
        Assumes that self.pos has been set correctly.
        """
        if self.pos is None:
            self.get_current_pos_vel()

        if self.field_values is None:
            self.field_values={}
        import string_unique
        fields_to_get = string_unique.unique(list(self.loop.fields_from_grid) + list(field_list))
        number_of_new_fields=0

        for field in fields_to_get:
            #initialize field values to -1.
            if field not in self.field_values:
                self.field_values[field] = np.zeros(self.ind.size)-1
                number_of_new_fields += 1

        if number_of_new_fields == 0:
            return

        verbose=False
        good_index_sort_np = np.array(copy.copy(self.ind)).astype('int64')
        gridlist = self.ds.index.grids[-1::-1]
        particle_fields = ['particle_position_x', 'particle_position_y', 'particle_position_z', 'particle_index', 'test_field']
        gotten = np.zeros(good_index_sort_np.size,dtype='int32') 
        for ngrid, grid in enumerate(gridlist):
            if verbose:
                print("grid %d / %d"%(ngrid, len(gridlist)))
            #grid i,j,k index selector.
            #mask between all particles and this grid.
            #mask_to_get_3 = np.zeros(good_index_sort_np.shape, dtype='int32')   #particles in this grid.  This is only for existence checking.
            #this gives the particles that live in _this grid._
            #found_any_g, mask_g = particle_ops.mask_particles(good_index_sort_np, grid['particle_index'].astype('int64'), mask_to_get_3)


            Nxyz = grid.ActiveDimensions
            pi =np.floor((self.pos[:,0] - grid.LeftEdge[0])/grid.dds[0]).astype('int32')
            pj =np.floor((self.pos[:,1] - grid.LeftEdge[1])/grid.dds[1]).astype('int32')
            pk =np.floor((self.pos[:,2] - grid.LeftEdge[2])/grid.dds[2]).astype('int32')
            pi=pi.v;pj=pj.v;pk=pk.v
            particle_selector = (pi >= 0 ) * ( pj >= 0) * ( pk >= 0 ) * ( pi < Nxyz[0])*(pj<Nxyz[1])*(pk<Nxyz[2])
            if particle_selector.sum() == 0:
                continue
            grid_selector = [pi,pj,pk]
            #this is giving a depreciation warning
            grid_to_particle = tuple([grid_selector[i][particle_selector] for i in [0,1,2]])
            subgrid_mask = grid.child_mask[grid_to_particle]
            particle_selector[particle_selector] = particle_selector[particle_selector]*subgrid_mask
            if particle_selector.max() > 0:  
                if verbose:
                    print("   work1")

                grid_to_particle_index = pk[particle_selector] + Nxyz[2]*(pj[particle_selector]+Nxyz[1]*pi[particle_selector])
                particle_sort_by_index = np.argsort(grid_to_particle_index)
                sorted_to_code_order = np.argsort(particle_sort_by_index)
                get_these_indices = grid_to_particle_index[particle_sort_by_index]
                particle_selector = particle_selector.astype('bool')
                for field in self.field_values:
                    if field[1] in particle_fields:
                        continue
                    #NOTE this mask can be move outside of the field loop.
                    if verbose:
                        print("   work1b")
                    arr = np.ascontiguousarray(grid[field].v).ravel()
                    values = arr[get_these_indices][sorted_to_code_order]
                    #values = grid[field][grid_to_particle]#grid_selector_tuple]
                    if verbose:
                        print("   work1c")
                    self.field_values[field][particle_selector] = values
                if verbose:
                    print("   work2")
                self.field_values[YT_particle_pos_x]=self.pos[:,0]
                self.field_values[YT_particle_pos_y]=self.pos[:,1]
                self.field_values[YT_particle_pos_z]=self.pos[:,2]
                self.field_values[YT_particle_index]=self.ind
                          
        #error check.
        if  (self.field_values[YT_cell_volume] < 0).any():
            NM= (self.field_values[YT_cell_volume] < 0).sum()
            print("ERROR: some particles (%d of them) not found.  This is problematic."%NM)
            pdb.set_trace()
        if verbose:
            print("   work3")

#Decorators frame_loop and core_loop takes functions and makes them loop over frame and core_id.
#The decorator takes care of loading the data, so your funciton only needs to use it.
#
#1.) write a function that does stuff as a method of the clump_looper
#2.) when you define it, put @looper_loop above the definition.
#3.) add it to the function
#4.) looper_loop shouldn't be modified
#See the 'full_plot' definition below,
#Also See decorator_test.py for more details.
def frame_loop(function):
    def wrapper(self,*args,**kwargs):
        for self.current_frame in self.frame_list:
            self.ds = self.load(self.current_frame)
            function(self,*args,**kwargs)
    return wrapper

def core_loop(function):
    def wrapper(looper,*args,**kwargs):
        for looper.current_frame in looper.frame_list:
            snapshot.ds = weakref.proxy(looper.load(frame=looper.current_frame))
            for looper.core_id in looper.core_list:
                this_snapshot = looper.make_snapshot(looper.current_frame,looper.core_id)
                if this_snapshot.R_centroid is None:
                    this_snapshot.get_all_properties()
                output = function(looper,this_snapshot,*args,**kwargs)
        return output #this is what the function will return
    return wrapper    #this return is from the decorator
# Decorators?
# A decorator looks like
# @some_decorator
# def some_function(arguments):
#    do_stuff()
# The decorator then does stuff to the way some_function behaves.
# So, if my decorator is do_five_times,
# def do_five_times(function):
#    def wrapper():
#        for i in range(5):
#            function()
# then I do
# @do_five_times
# def say_no():
#    print("no")
#
# say_no()
# will get called five times.  
# The decorator does
# say_no = do_five_times(say_no)
