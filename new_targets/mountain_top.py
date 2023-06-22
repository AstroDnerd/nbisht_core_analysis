from starter2 import *
import annotate_particles_4
from scipy.spatial import ConvexHull
from yt.data_objects.level_sets.api import *
from collections import defaultdict
import convex_hull_tools as CHT


class meta_locations():
    """Container object for data locations."""
    def __init__(self,this_simname,directory = None, frame = None, peak_fname=None, density_cut_fname=None):
        if directory is None:
            directory = dl.sims[this_simname]
        if frame is None:
            frame = dl.target_frames[this_simname]
        if peak_fname is None:
            peak_fname = dl.peak_list[this_simname]
        self.density_cuts={}
        if density_cut_fname is not None:
            dptr = h5py.File(density_cut_fname,'r')
            try:
                for n,peak_id in enumerate(dptr['peaks']):
                    density_cuts[peak_id]=dptr['density_cuts'][n]
            except:
                raise
            finally:
                dptr.close()
        ds_name="%s/DD%04d/data%04d"%(directory,frame,frame)
        print("OPENING %s and %s"%(ds_name, peak_fname))
        self.ds = yt.load(ds_name)
        fptr = h5py.File(peak_fname,'r')
        try:
            self.peaks = fptr['peaks'][()]
        except:
            raise
        finally:
            fptr.close()
        #peaks = nar([ np.array([0.89331055, 0.1159668 , 0.4440918 ])])

class top():
    """main mountain top object.  Data container, contains peak density and clumps."""
    def __init__(self, ds, location, top_to_bottom=3./4, rhomin=None, radius=1e-2, peak_id=-1):
        self.peak_id = peak_id
        self.radius = ds.arr(radius,'code_length')
        self.location =     ds.arr(location,'code_length')
        self.region = ds.sphere(self.location,self.radius)
        self.rhomax = get_density(self.location, self.region).v[0]
        self.rhomin = rhomin
        if self.rhomin is None:
             self.rhomin = self.rhomax**(top_to_bottom)
        self.clump = Clump(self.region,('gas','density'))
        self.clump.find_children( self.rhomin, self.rhomax)
        for leaf1 in self.clump.leaves:
            min_radius=leaf1['radius'].min()
            if min_radius < 1./2048:
                break
        self.leaf = leaf1

def check_overlap(leaf_storage):
    """Check for overlap between cores."""
    import convex_hull_tools as CHT
    from collections import defaultdict
    from scipy.spatial import ConvexHull
    overlap_list = defaultdict(list)
    leaves = np.array(list(leaf_storage.keys()))
    for nleaf,core_1 in enumerate(leaves):
        particles_1 = leaf_storage[core_1][1]['particle_position'].v
        if particles_1.shape[0] < 4:
            continue
        for core_2 in leaves[nleaf+1:]:
            #Is core_2 in core_1?
            print("check overlap on " , core_1, core_2)
            particles_2 = leaf_storage[core_2][1]['particle_position'].v
            hull = ConvexHull(particles_1)
            overlap = CHT.in_hull(particles_2,particles_1[hull.vertices,:])
            if overlap.any():
                print("====================== yes", core_1, core_2)
                overlap_list[core_1].append(core_2)
    return overlap_list

def is_point_in_hull(points2, test_main, tolerance=1./4096):
    """check of a point is in a convex hull of test_main"""
    X = test_main['x'].to('code_length')
    Y = test_main['y'].to('code_length')
    Z = test_main['z'].to('code_length')
    points1 = np.column_stack([X,Y,Z])

    try:
        hull1 = ConvexHull( points1 )
        in_hull_1 = CHT.in_hull( [points2], points1[hull1.vertices,:])
    except:
        print("Cannot find hull.  Assume exclusion.")
        in_hull_1 = False
    return in_hull_1

def leaf_with_center(leaf_list):
    """This has a hard coded minimum point size.  For finding child leaves, maybe not the best method."""
    for leaf in leaf_list:
        if leaf['radius'].min()<1/4096:
            return leaf
def get_density(peak, test_main):
    """Find the density at a location within a yt region"""
    X = test_main[YT_x].to('code_length')
    Y = test_main[YT_y].to('code_length')
    Z = test_main[YT_z].to('code_length')
    tolerance = 0.5/2048
    mask = ( np.abs(X - peak[0])<tolerance)*( np.abs(Y-peak[1])<tolerance)*(np.abs(Z-peak[2])<tolerance)
    if mask.sum() != 1:
        pdb.set_trace()
    return test_main['density'][mask]

def read_peaks(this_simname):
    """read the peak hdf5 file"""
    peak_fname = dl.peak_list[this_simname]
    fptr = h5py.File(peak_fname,'r')
    try:
        peaks = fptr['peaks'][()]
    except:
        raise
    finally:
        fptr.close()
    return peaks

class target_info():
    """Target object for the tracker."""
    def __init__(self,  peak_id = None, peak_density=None, min_density=None,
                 peak_location = None, nzones=None,
                 particle_index = None, h5ptr=None):
        if h5ptr == None:
            self.peak_id = peak_id
            self.peak_density=  peak_density
            self.min_density= min_density
            self.peak_location = peak_location
            self.nzones=nzones
            self.particle_index = particle_index
        else:
            self.peak_id        =h5ptr["peak_id"][()]
            self.peak_density   =h5ptr["peak_density"][()]
            self.min_density    =h5ptr["min_density"][()]
            self.peak_location  =h5ptr["peak_location"][()]
            self.nzones         =h5ptr["nzones"][()]
            self.particle_index =h5ptr["particle_index"][()]

    def write(self,h5ptr):
        group_name = "peak_p%04d"%self.peak_id
        group = h5ptr.create_group(group_name)
        group.create_dataset("peak_id",       data = self.peak_id)
        group.create_dataset("peak_density",  data = self.peak_density)
        group.create_dataset("min_density",   data = self.min_density)
        group.create_dataset("peak_location", data = self.peak_location)
        group.create_dataset("nzones",        data = self.nzones)
        group.create_dataset("particle_index",data = self.particle_index)



def cut_mountain_top(this_simname, target_fname, density_cut_fname=None, directory = None, frame = None, peak_fname=None, work_radius=1e-2, do_projections=False, top_to_bottom=3./4, kludge={},
                    verify=None, leaf_storage=None,radius=1e-2, cut_override={}, radius_dict={}):

    """Loop over peaks.  Make mountain tops."""
    meta = meta_locations(this_simname,directory=directory,frame=frame,peak_fname=peak_fname,density_cut_fname=density_cut_fname)
    meta.density_cuts.update(cut_override)

    ds = meta.ds
    h5ptr = h5py.File(target_fname,'w')
    for peak_id,center_o in enumerate(meta.peaks):
        if "peak_id" in kludge:
            if peak_id not in kludge['peak_id']:
                continue


        radius = ds.arr(radius_dict.get(peak_id,1e-2),'code_length')
        top1 = top(ds,meta.peaks[peak_id], radius=radius,top_to_bottom=top_to_bottom, rhomin=meta.density_cuts.get(peak_id,None), peak_id=peak_id)
        this_target_info = target_info( peak_id = peak_id, peak_density=top1.rhomax, min_density=top1.rhomin,
                                       peak_location = center_o, nzones=top1.leaf['radius'].size,
                                       particle_index = np.sort(top1.leaf['particle_index']))
        FAIL = ""
        write_dataset=True
        if verify is not None and not verify(top1):
            print("CUTTING DATSET")
            write_dataset=False
            FAIL = "FAIL "


        if leaf_storage is not None and type(leaf_storage) == dict and write_dataset:
            leaf_storage[peak_id] = top1,top1.leaf


        if write_dataset:
            print("WRITE")
            try:
                this_target_info.write(h5ptr)
            except:
                #it's really irritating if you accidentally leave an hdf5 file open.
                h5ptr.close()
                raise
        if do_projections:
            proj = ds.proj('density',0,center=top1.location,data_source=top1.region)
            pw = proj.to_pw()
            pw.set_cmap('density','Greys')

            #pw.annotate_clumps([master_clump]+master_clump.leaves)
            if top1.leaf['particle_index'].size > 10:
                p_size = 1
            else:
                p_size = 7
            pw.annotate_these_particles4(1.0,col='r',positions= top1.leaf['particle_position'], p_size=p_size)
            pw.zoom(0.5/radius.v)
            pw.set_axes_unit('code_length')

            pw.annotate_clumps([top1.leaf])
            pw.annotate_title("%s Peak id %d N particles %d in leaf %d rho %0.2e"%(FAIL,peak_id,top1.region['particle_index'].size, top1.leaf['particle_index'].size, top1.rhomax ))
            pw.save('plots_to_sort/%s_peak_p%04d_34'%(this_simname,peak_id))
            #            proj = ds.proj(field,ax,center=center, data_source = sph) 
    h5ptr.close()


def split_peaks(this_simname, peak_1, peak_2,directory = None, frame = None, peak_fname=None, density_cut_fname=None,
                work_radius=1e-2, do_projections=False, top_to_bottom=3./4, kludge={},
                recheck=False,   verify=None, leaf_storage=None, radius_dict={}, auto_double=True,contour_step=1.5):
    meta = meta_locations(this_simname,directory=directory,frame=frame,peak_fname=peak_fname,density_cut_fname=density_cut_fname)
    ds = meta.ds
    radius1 = ds.arr(radius_dict.get(peak_1,1e-2),'code_length')
    radius2 = ds.arr(radius_dict.get(peak_2,1e-2),'code_length')
    print("RADIUS", radius1, radius2)
    output={}

    #NOTE need to pass in alternate contours.
    top1 = top(ds,meta.peaks[peak_1],radius=radius1,top_to_bottom=top_to_bottom,peak_id=peak_1)
    top2 = top(ds,meta.peaks[peak_2],radius=radius2,top_to_bottom=top_to_bottom,peak_id=peak_2)

    #low density pairs need to be bigger
    if top1.rhomax < 3000 and top2.rhomax < 3000 and auto_double:
        print("BOTH LOW DENSITY USING BIGGER SPHERE")
        radius1 = ds.arr(radius_dict.get(peak_1,2e-2),'code_length')
        radius2 = ds.arr(radius_dict.get(peak_2,2e-2),'code_length')
        top1 = top(ds,meta.peaks[peak_1],radius=radius1,top_to_bottom=top_to_bottom,peak_id=peak_1)
        top2 = top(ds,meta.peaks[peak_2],radius=radius2,top_to_bottom=top_to_bottom,peak_id=peak_2)
    


    which_to_contour=np.argmin([top1.rhomax,top2.rhomax])
    which_other     =np.argmax([top1.rhomax,top2.rhomax])
    tops = [top1,top2]
    top_to_contour = tops.pop(which_to_contour)
    top_other = tops[0]
    peak_to_contour = [peak_1,peak_2][which_to_contour]
    peak_other      = [peak_1,peak_2][which_other]

    leaf_storage['orig']={}
    leaf_storage['orig'][peak_1]=ds,top1.leaf
    leaf_storage['orig'][peak_2]=ds,top2.leaf

    points1 = top1.leaf['particle_position'].v
    points2 = top2.leaf['particle_position'].v
    hull1 = ConvexHull( points1)
    try_anyway=True
    try:
        in_hull_1 = CHT.in_hull( points2, points1)#points1[hull1.vertices,:])
        n_overlap = in_hull_1.sum()
    except:
        try_anyway = True
        n_overlap=0
    do_child_plots=True
    if n_overlap == 0 and not try_anyway:
        print("ERROR no overlap, %d %d"%(peak_1,peak_2))
        output[peak_1]=top1.rhomin
        output[peak_2]=top2.rhomin
        do_child_plots = False


    else:
        min_density = top_to_contour.rhomin
        temp_clump = Clump(top_to_contour.region, ('gas','density'), clump_info="")
        temp_clump.add_validator("min_cells", 8)
        find_clumps(temp_clump, top_to_contour.rhomin, top_to_contour.rhomax, contour_step)

        #Find the set of clumps that contain the center of the region 
        this_stack = []
        loop_clump = temp_clump
        while loop_clump.children is not None and len(loop_clump.children):
            central_leaf =  leaf_with_center( loop_clump.children)
            if central_leaf is not None:
                this_stack.append( central_leaf )
            loop_clump = this_stack[-1]
        new_clump=None
        #Find the first clump that rejects the other point
        all_excluded=[]
        point_to_exclude = top_other.location
        for clump in this_stack:
            c2_in_new_clump= is_point_in_hull( point_to_exclude, clump) 
            if not c2_in_new_clump:
                all_excluded.append( clump )
        if len(all_excluded):
            new_clump = all_excluded[0]
        else:
            print("Cannot find excluded contour.  Please help me.")
            new_clump = None
            do_child_plots = False
        if new_clump is None:
            print("Cannot find a new clump.  Reject one of them.")
            output[peak_1]=-1
            output[peak_2]=-1
            do_child_plots=False

        else:
            if top_to_contour.rhomin > top_other.rhomax:
                output[peak_to_contour ] =max([top_to_contour.rhomin, new_clump.contour_level])
                output[peak_other ]      =max([top_other.rhomin,      top_other.rhomin])
            else:                         
                output[peak_to_contour ] =max([top_to_contour.rhomin, new_clump.contour_level])
                output[peak_other ]      =max([top_other.rhomin,      new_clump.contour_level])





        leaf_storage['new_clump']=new_clump
        leaf_storage['excluded'] = all_excluded
        leaf_storage['temp_clump'] = temp_clump
        leaf_storage['stuff']={'top1':top1,'top2':top2}

    if recheck:
        top3 = top(ds,meta.peaks[peak_1], radius1, rhomin=output[peak_1],peak_id=peak_1)
        top4 = top(ds,meta.peaks[peak_2], radius2, rhomin=output[peak_2],peak_id=peak_2)
        leaf_storage['new']={}
        leaf_storage['new'][peak_1]=ds,top3.leaf
        leaf_storage['new'][peak_2]=ds,top4.leaf


    verbose=True
    if verbose:
        print("ZZZ Core 1:  rhomax %0.2e rhomin %0.2enew %0.2e"%( top1.rhomax, top1.rhomin,output[peak_1]))
        print("ZZZ Core 2:  rhomax %0.2e rhomin %0.2enew %0.2e"%( top2.rhomax, top2.rhomin,output[peak_2]))
    if do_projections:
        proj = ds.proj('density',0,center=top1.location,data_source=top1.region)
        pw = proj.to_pw()
        pw.set_cmap('density','Greys')
        pw.set_axes_unit('code_length')
        pw.zoom(0.5/radius1.v)

        #pw.annotate_clumps([master_clump]+master_clump.leaves)
        p_size=1
        pw.annotate_these_particles2(1.0,col='r',positions= top1.leaf['particle_position'], p_size=p_size)
        pw.annotate_these_particles2(1.0,col='g',positions= top2.leaf['particle_position'], p_size=p_size)
        pw.annotate_clumps([leaf_storage['orig'][peak_1][1], leaf_storage['orig'][peak_2][1]])
        pw.annotate_text(top1.location,'*')
        pw.annotate_text(top2.location,'x')
        pw.save('plots_to_sort/%s_peak_split_orig_c%04d_c%04d'%(this_simname,peak_1,peak_2))

        #pw.annotate_clear()
        #pw.annotate_clumps(temp_clump.leaves)
        #pw.annotate_text(top1.location,'*')
        #pw.annotate_text(top2.location,'x')
        #print(c1,c2)
        #pw.save('plots_to_sort/%s_peak_split_temp_leaves_c%04d_c%04d'%(this_simname,peak_1,peak_2))

        pw.clear_annotations()
        pw.annotate_clumps(this_stack[::-1])
        pw.annotate_text(top1.location,'*')
        pw.annotate_text(top2.location,'x')
        pw.save('plots_to_sort/%s_peak_split_corestack_c%04d_c%04d'%(this_simname,peak_1,peak_2))

        if do_child_plots:
            pw.clear_annotations()
            pw.annotate_clumps(all_excluded[::-1])
            pw.annotate_text(top1.location,'*')
            pw.annotate_text(top2.location,'x')
            pw.save('plots_to_sort/%s_peak_split_all_excluded_c%04d_c%04d'%(this_simname,peak_1,peak_2))
            pw.clear_annotations()
            pw.annotate_clumps([new_clump])
            pw.annotate_text(top1.location,'*')
            pw.annotate_text(top2.location,'x')
            pw.save('plots_to_sort/%s_peak_split_newclump_c%04d_c%04d'%(this_simname,peak_1,peak_2))
            if recheck:
                pw.clear_annotations()
                pw.annotate_clumps([top3.leaf,top4.leaf])
                pw.annotate_text(top1.location,'*')
                pw.annotate_text(top2.location,'x')
                pw.save('plots_to_sort/%s_peak_split_leaf3leaf4_c%04d_c%04d'%(this_simname,peak_1,peak_2))
            #            proj = ds.proj(field,ax,center=center, data_source = sph) 
    return output

def split_all(this_simname, overlap_all, new_thresholds, all_peaks, skip_list=[], peak_list=None, do_projections=False, radius_dict={},
              auto_double=True,contour_step=1.5):
    if peak_list is None:
        peak_list = sorted( list(overlap_all.keys()))
    for peak_id in peak_list:
        if peak_id in skip_list:
            continue
        print("=============================================", peak_id)
        for other_peak in overlap_all[peak_id]:
            if other_peak in skip_list:
                continue
            if peak_id in all_peaks:
                if other_peak in all_peaks[peak_id]:
                    continue
            print("DO COMPARE", peak_id, other_peak)
            radius = radius_dict.get(peak_id,1e-2)
            more_leaf_storage={}
            split_density=split_peaks(this_simname, peak_id,other_peak ,do_projections=do_projections,
                                      leaf_storage=more_leaf_storage,radius_dict=radius_dict, auto_double=auto_double, 
                                      contour_step=contour_step)
            all_peaks[peak_id][other_peak]=split_density[peak_id]
            for ppp in split_density:
                this_density = split_density[ppp]
                new_thresholds[ppp] = max([this_density,new_thresholds.get(ppp,this_density)])

#
# old tools, keep around until we're sure we're done.
#
def plot_leaves(this_simname, target_fname, density_cut_fname=None, directory = None, frame = None, peak_fname=None, work_radius=1e-2, do_projections=False, top_to_bottom=3./4, kludge={},
                    verify=None, leaf_storage=None):
    if directory is None:
        directory = dl.sims[this_simname]
    if frame is None:
        frame = dl.target_frames[this_simname]
    if peak_fname is None:
        peak_fname = dl.peak_list[this_simname]
    density_cuts={}
    if density_cut_fname is not None:
        dptr = h5py.File(density_cut_fname,'r')
        try:
            for n,peak_id in enumerate(dptr['peaks']):
                density_cuts[peak_id]=dptr['density_cuts'][n]
        except:
            raise
        finally:
            dptr.close()
    ds_name="%s/DD%04d/data%04d"%(directory,frame,frame)
    print("OPENING %s and %s"%(ds_name, peak_fname))
    ds = yt.load(ds_name)

def get_peak_densities(this_simname, target_fname, density_cut_fname=None, directory = None, frame = None, peak_fname=None, work_radius=1e-2, do_projections=False, top_to_bottom=3./4, kludge={},
                    verify=None, leaf_storage=None):
    if directory is None:
        directory = dl.sims[this_simname]
    if frame is None:
        frame = dl.target_frames[this_simname]
    if peak_fname is None:
        peak_fname = dl.peak_list[this_simname]
    density_cuts={}
    if density_cut_fname is not None:
        dptr = h5py.File(density_cut_fname,'r')
        try:
            for n,peak_id in enumerate(dptr['peaks']):
                density_cuts[peak_id]=dptr['density_cuts'][n]
        except:
            raise
        finally:
            dptr.close()
    ds_name="%s/DD%04d/data%04d"%(directory,frame,frame)
    print("OPENING %s and %s"%(ds_name, peak_fname))
    ds = yt.load(ds_name)
    fptr = h5py.File(peak_fname,'r')
    try:
        peaks = fptr['peaks'][()]
    except:
        raise
    finally:
        fptr.close()
    #peaks = nar([ np.array([0.89331055, 0.1159668 , 0.4440918 ])])
    radius = ds.arr(1e-2,'code_length')
    peak_densities=[]
    for peak_id,center_o in enumerate(peaks):
        if "peak_id" in kludge:
            if peak_id not in kludge['peak_id']:
                continue

        center = ds.arr(center_o, 'code_length')
        test_main = ds.sphere(center, radius)
        peak_density = get_density(center,test_main)
        peak_densities.append(peak_density)
        print(" peak %d %0.2e"%(peak_id,peak_density))
    return peak_densities

def dump_cutoffs(fname,split_dict):
    peaks=[]
    density_cuts=[]
    if os.path.exists(fname):
        fptr = h5py.File(fname,'r')
        peaks = list(fptr['peaks'][()])
        density_cuts = list(fptr['density_cuts'][()])
        fptr.close()
    for peak_id in split_dict:
        if peak_id in peaks:
            ind = peaks.index(peak_id)
            density_cuts[ind] = max( density_cuts[ind], split_dict[peak_id])
        else:
            peaks.append(peak_id) 
            density_cuts.append(split_dict[peak_id])
    fptr=h5py.File(fname,'w')
    fptr['peaks']=peaks
    fptr['density_cuts']=density_cuts
    fptr.close()
def recursive_print(children):
    for child in children:
        print(child.contour_level)
        if len(child.children):
            recursive_print(child.children)

class target_tree():
    """Unused."""
    def __init__(self,ds=None,fname=None):
        self.targets=None
        self.target_indices={}
        self.leaves={}
        self.ds=ds
        if fname is not None:
            self.read_targets(fname)
    def read_targets(self,target_fname):
        self.all_particles=np.array([])
        self.core_ids_by_particle=np.array([])
        if self.targets is None:
            self.targets = {}
        h5ptr = h5py.File(target_fname,'r')
        for group_name in h5ptr:
            core_id = h5ptr[group_name]['peak_id'][()]
            #read from disk
            self.targets[core_id] = target_info(h5ptr=h5ptr[group_name])
            self.target_indices[core_id] = self.targets[core_id].particle_index

            #check for uniqueness
            self.all_particles=np.append(self.all_particles, self.target_indices[core_id] )
            if np.unique(self.all_particles).size - self.all_particles.size:
                print("FATAL ERROR: repeated particle, ", core_id)
                pdb.set_trace()
            #I might need this later when I revamp get_tracks again.
            these_core_ids=[core_id]*self.target_indices[core_id].size
            self.core_ids_by_particle=np.append(self.core_ids_by_particle, these_core_ids)

        h5ptr.close()
        self.core_list = np.sort( np.unique( list(self.target_indices.keys())))
    def make_leaves_check(self):
        if self.leaves is None:
            self.leaves={}
        for core_id in self.targets:
            T = self.targets[core_id]
            sph = self.ds.sphere( T.peak_location, 1e-2)
            clump = Clump(sph,('gas','density'))
            clump.find_children( T.min_density, T.peak_density)
            for leaf1 in clump.children:
                min_radius=leaf1['radius'].min()
                if min_radius < 1./2048:
                    break
            self.leaves[core_id] = leaf1

            self.leaves
