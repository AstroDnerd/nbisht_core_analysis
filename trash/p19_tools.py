import yt
import matplotlib.pyplot as plt
ef=execfile
from yt.visualization.plot_modifications import *
import pyximport; pyximport.install()
import particle_ops
import particle_grid_mask
import fPickle
import time
import numpy as na
import os
import pdb
import copy
dbg = 0
def get_leaf_indices(ds,c_min=None,c_max=None,step=100,pickle_name=None, subset=None):
    """get all the leaf indices for peaks in *ds*.
    If *pickle_name* is supplied, load from that, or if it doesn't exist, save to that.
    *subset*, if supplied, will restrict the indices returned.
    """
    if not os.path.exists(pickle_name):
        print("WARNING: clump finding not tested.")
        ad = ds.all_data()
        master_clump = Clump(ad,('gas','density'))
        master_clump.add_validator("min_cells", 8)
        c_min = ad["gas", "density"].min()
        c_max = ad["gas", "density"].max()
        step = 100
        find_clumps(master_clump, c_min, c_max, step)
    # Write a text file of only the leaf nodes.
        write_clumps(master_clump,0, "%s_clumps.txt" % ds)

        leaf_clumps = get_lowest_clumps(master_clump)
    

        for i in range(len(leaf_clumps)):
            den_max.append(leaf_clumps[i][('gas','density')].max())
            x_max.append(leaf_clumps[i]['x'][np.where(leaf_clumps[i]['gas','density']==den_max[i])])
            y_max.append(leaf_clumps[i]['y'][np.where(leaf_clumps[i]['gas','density']==den_max[i])])
            z_max.append(leaf_clumps[i]['z'][np.where(leaf_clumps[i]['gas','density']==den_max[i])])
  
            a= float(x_max[i])
            b= float(y_max[i])
            c= float(z_max[i])
            this_peak = ds.arr([a,b,c],'code_length')
            peak_list.append(loc)
        fPickle.dump(peak_list, pickle_name)
    else:
        peak_list = fPickle.load(pickle_name)
    
    if subset is None:
        subset = range(len(peak_list))

    leaf_indices={}

    min_dx = ds.index.get_smallest_dx()
    for clump in subset:
        this_clump = peak_list[clump]
        if not hasattr(this_clump,'units'):
            this_clump = ds.arr(this_clump,'code_length')

        region = ds.region(center=this_clump,left_edge = this_clump-1.5*min_dx, right_edge = this_clump+1.5*min_dx)
        indices = region['particle_index'].astype('int64')
        leaf_indices[clump]=indices
    return leaf_indices
    #for nc,indices in enumerate(leaf_indices):
     #   pw_full.annotate_select_particles(1.0, col='r', indices=indices)
   # pw_full.save(fname)
def shift_particles(ds, position,shiftRight = False):
    """Shifts a periodically separated clump by the domain width.
    Looks for gaps in the positions larger than max('dx'), shifts one group
    to the right (left if shiftRight=Flase) to be spatially contiguous."""
    #max_dx may not be computed in the most efficient way.
    DomainLeft = ds.domain_left_edge
    DomainRight =  ds.domain_right_edge
    DomainWidth = DomainRight - DomainLeft
    shift = np.zeros(3)
    shifted=copy.copy(position)
    for i,axis in enumerate(['x','y','z']):

        dx = 'd'+axis
        nique = np.unique(shifted[:,i])
        nique.sort()
        max_dx = ds.index.grids[0].dds.max().in_units('code_length')
        min_dx = ds.index.get_smallest_dx()

        #has to be close to the edges, or 'Periodic Wrap' isn't the problem.
        if np.abs(nique.max() - DomainRight[i]) > 3*max_dx:
            continue
        if np.abs(nique.min() - DomainLeft[i]) > 3*max_dx:
            continue
        delta_x = nique[1:] - nique[0:-1]
        max_delta_x = delta_x.max()
        break_index = np.where(delta_x == max_delta_x)
        if max_delta_x > max_dx:
           if shiftRight:
              break_x = nique[break_index[0]]
              all_to_shift = np.where( shifted[:,i] <= break_x + min_dx )[0]
              shifted[:,i][all_to_shift] += DomainWidth[i]
              shift[i] = DomainWidth[i]
           else:
              break_x = nique[break_index[0]+1]
              all_to_shift = np.where( shifted[:,i] >= break_x - min_dx )[0]
              shifted[:,i][all_to_shift] -= DomainWidth[i]
              shift[i] = -DomainWidth[i]
                
    return shifted
def roundup(x):
    return int(math.ceil(x/10.0))*10
def powerline(this_plt,x1,x2,y1,power,**kwargs):
    """Plot a powerlaw on the current plot in matplot lib instance *plt*.
    Starts at *x1*, *y1* and runs to *x2*, and the y2 position dictated by the exponent.*log* determines log or linear plot."""
    x = [x1,x2]
    yf = x2**power*x1**(-power)*y1
    y = [y1,yf]
    this_plt.plot(x,y,**kwargs)

def deposit_target_particles_1(field, data):
    #pdb.set_trace()
    indices_late = data.get_field_parameter('indices_late')
    blank = np.zeros(data.ActiveDimensions, dtype='float64')
    if not hasattr(indices_late,'sum'):
        return blank
    if data.NumberOfParticles == 0: return blank
    if dbg > 0:
        if hasattr(indices_late,'sum'):
            print "indices inside", indices_late.sum()
        else:
            print "indices inside", indices_late
            return blank
    #int 32 because it's just a flag.
    #mask_to_get = np.zeros(indices_late.shape, dtype='int32')
    mask_to_get = data.get_field_parameter('mask_to_get')
    t0 = data.get_field_parameter('timer')
    my_indices = data['particle_index'].astype('int64')
    #print "  min thing"
    #mask_to_get[ indices_late < my_indices.min()] = 1
    #mask_to_get[ indices_late > my_indices.max()] = 1
    #print "  left", mask_to_get.size - mask_to_get.sum()
    #print "  search"
    found_any, mask = particle_ops.mask_particles(
        indices_late, my_indices, mask_to_get)
    data.set_field_parameter('mask_to_get',mask_to_get)
    pos_to_get = data['particle_position'][mask == 1]
    d = data.deposit(pos_to_get, method = "count")
    d = data.ds.arr(d, input_units = "cm**-3")
    t1 = time.time() 
    #print "  DT", t1-t0[-1]
    data.set_field_parameter('timer',t0+[t1])
    return data.apply_units(d, field.units)

yt.add_field(      ("deposit","deposit_target_particles"),
#yt.add_field(      'dp1',
         function = deposit_target_particles_1,
         validators = [yt.ValidateSpatial(), 
                       yt.ValidateParameter('indices_late'), 
                       yt.ValidateParameter('mask_to_get'), 
                       yt.ValidateGridType()],
         display_name = "target_particles")

class SelectParticleCallback(PlotCallback):
    _type_name = "select_particles"
    region = None
    _descriptor = None
    def __init__(self, width, p_size=1.0, col='k', marker='o', stride=1.0,
                 ptype='all', stars_only=False, dm_only=False,
                 minimum_mass=None,bool=None, indices=None, xyz=None):
        """
        Adds particle positions, based on a thick slab along *axis* with a
        *width* along the line of sight.  *p_size* controls the number of
        pixels per particle, and *col* governs the color.  *ptype* will
        restrict plotted particles to only those that are of a given type.
        *minimum_mass* will require that the particles be of a given mass,
        calculated via ParticleMassMsun, to be plotted.
        """
        PlotCallback.__init__(self)
        self.width = width
        self.p_size = p_size
        self.color = col
        self.marker = marker
        self.stride = stride
        self.ptype = ptype
        self.stars_only = stars_only
        self.dm_only = dm_only
        self.minimum_mass = minimum_mass
        self.bool = bool
        self.indices = indices
        self.xyz = xyz

    def __call__(self, plot):
        data = plot.data
        if iterable(self.width):
            self.width = np.float64(plot.data.ds.quan(self.width[0], self.width[1]))
        # we construct a recantangular prism
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        if self.bool is None:
            reg = self._get_region((x0,x1), (y0,y1), plot.data.axis, data)
        else:
            reg = self.bool
        ax = data.axis


        xax = plot.data.ds.coordinates.x_axis[ax]
        yax = plot.data.ds.coordinates.y_axis[ax]
        axis_names = plot.data.ds.coordinates.axis_name
        field_x = "particle_position_%s" % axis_names[xax]
        field_y = "particle_position_%s" % axis_names[yax]
        pt = self.ptype
        shifted_field_x = reg[pt, field_x]
        shifted_field_y = reg[pt, field_y]
        shifted_field = [shifted_field_x, shifted_field_y]
        lower_lim = data.ds.arr([xx0,yy0],'code_length')
        upper_lim = data.ds.arr([xx1,yy1],'code_length')
        dx_min = data.ds.index.get_smallest_dx()*0.1

        for ndim,dim in enumerate([xax,yax]):
            if data.ds.periodicity[dim] and True:
                print("SHIFT", ndim, dim)
                #If the plot is to the right of the domain
                delta = upper_lim[ndim] - data.ds.domain_right_edge[dim]
                print("DELTA up", dim, delta, upper_lim[ndim])
                if delta > dx_min:
                    delta = delta + data.ds.domain_left_edge[dim]
                    particles_to_shift = np.where(shifted_field[ndim] < delta)
                    shifted_field[ndim][particles_to_shift] += data.ds.domain_width[dim]
                #If the plot is to the left
                delta = data.ds.domain_left_edge[dim] - lower_lim[ndim]
                print("DELTA down", dim, delta, lower_lim[ndim])
                if delta > dx_min:
                    delta = data.ds.domain_right_edge[dim] - delta
                    particles_to_shift = np.where(shifted_field[ndim] > delta)
                    shifted_field[ndim][particles_to_shift] -= data.ds.domain_width[dim]



        gg = ( ( shifted_field_x >= x0 ) & ( shifted_field_x <= x1 )
           &   ( shifted_field_y >= y0 ) & ( shifted_field_y <= y1 ) )
        print("CHECK", x0, x1, y0, y1, "npart", reg['particle_index'].size,"n(gg)", gg.sum())
        print("CHECK", shifted_field_x.min(), shifted_field_x.max(), shifted_field_y.min(),shifted_field_y.max())
        #import pdb
        #pdb.set_trace()

        if self.indices != None:
            print "particle callback warning: this might get expensive"
            mask_to_get = na.zeros(self.indices.shape, dtype='int32')
            found_any, mask = particle_ops.mask_particles(
                self.indices.astype('int64'), reg['particle_index'].astype('int64'), mask_to_get)
            gg = ( gg & (mask == 1) )

        print "nparticles do it.", mask_to_get.sum()
        if False:
            if self.ptype is not None:
                gg &= (reg["particle_type"] == self.ptype)
                if gg.sum() == 0: return
            if self.stars_only:
                gg &= (reg["creation_time"] > 0.0)
                if gg.sum() == 0: return
            if self.dm_only:
                gg &= (reg["creation_time"] <= 0.0)
                if gg.sum() == 0: return
            if self.minimum_mass is not None:
                gg &= (reg["ParticleMassMsun"] >= self.minimum_mass)
                if gg.sum() == 0: return
        #plot._axes.hold(True)
        px, py = self.convert_to_plot(plot,
                    [np.array(shifted_field_x[gg][::self.stride]),
                     np.array(shifted_field_y[gg][::self.stride])])
        plot._axes.scatter(px, py, edgecolors='None', marker=self.marker,
                           s=self.p_size, c=self.color)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        #plot._axes.hold(False)

    def _get_region(self, xlim, ylim, axis, data):
        LE, RE = [None]*3, [None]*3
        ds=data.ds
        xax = ds.coordinates.x_axis[axis]
        yax = ds.coordinates.y_axis[axis]
        zax = axis
        LE[xax], RE[xax] = xlim
        LE[yax], RE[yax] = ylim
        LE[zax] = data.center[zax].ndarray_view()- self.width*0.5
        RE[zax] = data.center[zax].ndarray_view()+ self.width*0.5
        if self.region is not None \
            and na.all(self.region.left_edge <= LE) \
            and na.all(self.region.right_edge >= RE):
            return self.region
        self.region = data.ds.region( data.center, LE, RE)
        return self.region

def get_indices(these_pids, ad):
    these_pids=these_pids.astype('int64').v
    mask_to_get=np.zeros(these_pids.shape,dtype='int32')
    my_indices = ad['particle_index'].astype('int64')
    found_any, mask = particle_ops.mask_particles(these_pids,my_indices,mask_to_get)
    return found_any, mask
    #good_index = ad['particle_index'][mask == 1]
    #good_index_sort = sorted(good_index,key = int)
    #good_index_sort_np = np.array(good_index_sort).astype('int64')
    #pos = ad['particle_position'][mask == 1]
    #vel = ad['particle_velocity'][mask == 1]

def get_particle_values_from_grid(ds,field_values,leaf_indices, particle_pos, 
                                  individual_particle_tracks=False,
                                  pdict = {}):

    good_index_sort_np = np.array(copy.copy(leaf_indices)).astype('int64')
    for grid in ds.index.grids[-1::-1]:
        grid_selector = np.zeros([3,good_index_sort_np.size],dtype='int32') #grid i,j,k index selector.
        particle_selector = np.zeros(good_index_sort_np.size,dtype='int32') #mask between all particles and this grid.
        mask_to_get_3 = np.zeros(good_index_sort_np.shape, dtype='int32')   #particles in this grid.  This is only for existence checking.
        #this gives the particles that live in _this grid._
        found_any_g, mask_g = particle_ops.mask_particles(good_index_sort_np, grid['particle_index'].astype('int64'), mask_to_get_3)
        if found_any_g:
            particle_grid_mask.particle_grid_mask_go_i2(particle_pos[:,0],particle_pos[:,1],particle_pos[:,2], 
                                                        grid.LeftEdge, grid.dds, 
                                                        grid.ActiveDimensions,grid.child_mask, 
                                                        grid_selector,
                                                        particle_selector)

            particle_selector = particle_selector==1
            if particle_selector.sum() == 0:
                continue
            for field in field_values.keys():
                values = grid[field][[grid_selector[i][particle_selector] for i in [0,1,2]]]
                field_values[field][particle_selector] = values
            if individual_particle_tracks:
                particles_in_grids = these_pids[particle_selector]
                for p in particles_in_grids:
                   pdict[p]=pdict.get(p,{})
                   for k in ['times','cycles']:
                       if not pdict[p].has_key(k): pdict[p][k]=[]
                   pdict[p]['times'].append(times[-1])
                   pdict[p]['cycles'].append(cycles[-1])
                for field in ['density']:
                   values = grid[field][[grid_selector[i][particle_selector] for i in [0,1,2]]]
                   for n,p in enumerate(particles_in_grids):
                       pdict[p][field] = pdict[p].get(field,[])
                       pdict[p][field].append(values[n])
