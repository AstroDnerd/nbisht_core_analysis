import yt
from tools.yt_names import *
import matplotlib.pyplot as plt
from yt.visualization.plot_modifications import *
import pyximport; pyximport.install()
from scipy.spatial import ConvexHull
import h5py
import time
import numpy as na
import os
import pdb
import copy
dbg = 0
from yt.data_objects.level_sets.api import *

def check_particles(ds):
    bad_index=[]
    for grid in ds.index.grids:
        pos = grid['particle_position']
        for i in [0,1,2]:
            check=pos[:,i] > grid.RightEdge[i]

            check=np.logical_or(check,pos[:,i] < grid.LeftEdge[i])
            if check.any():
                locations=np.where(check)
                particles= grid['particle_index'][check].v
                bad_index+= list(particles)
    return bad_index
#g=check_particles(this_looper.ds_list[120])
def get_leaf_clumps(ds,c_min=None,c_max=None,step=None, small_test=False, density_placeholder = 'density'):
    """get all the leaf indices for peaks in *ds*.
    """
    if small_test:
        #center = ds.arr([0.07104492, 0.05688477, 0.1862793 ],'code_length')
        peak,center=ds.find_max(density_placeholder)
        ad = ds.sphere(center,0.1)
    else:
        ad = ds.all_data()
    #ad  = ds.sphere([0.52075195, 0.74682617, 0.01196289], 0.1)
    master_clump = Clump(ad,('gas',density_placeholder))
    master_clump.add_validator("min_cells", 8)
    if c_min is None:
        c_min = 10 #ad["gas", "density"].min()
    #c_max = 534069645. # ad["gas", "density"].max()
    if c_max is None:
        c_max = ad["gas", density_placeholder].max()
    if step is None:
        step = 100
    find_clumps(master_clump, c_min, c_max, step)
    return master_clump

from yt.data_objects.level_sets.clump_tools import return_bottom_clumps
def get_peak_indices(master_clump,ds,h5_name="file.h5", density_placeholder = 'density'):

    leaf_clumps = return_bottom_clumps(master_clump)

    peak_list=[]
    den_max=[]
    x_max=[]
    y_max=[]
    z_max=[]
    for i in range(len(leaf_clumps)):
        den_max.append(leaf_clumps[i][('gas',density_placeholder)].max())
        max_loc = np.where(leaf_clumps[i]['gas',density_placeholder]==den_max[i])
        a = leaf_clumps[i][YT_x][max_loc][0]
        b = leaf_clumps[i][YT_y][max_loc][0]
        c = leaf_clumps[i][YT_z][max_loc][0]

        this_peak = ds.arr([a,b,c],'code_length')
        peak_list.append(this_peak)
    #fPickle.dump(peak_list, pickle_name)
    fptr=h5py.File(h5_name,'w')
    fptr.create_dataset('peaks',data=np.array(peak_list))
    fptr.close()
def get_leaf_indices(ds,c_min=None,c_max=None,step=100,h5_name="NEW_PEAK_FILE.h5",pickle_name=None, 
                     subset=None, peak_radius=1.5,bad_particle_list=None, small_test=False):

    if 1:
        fptr = h5py.File(h5_name,'r')
        peak_list = fptr['peaks'][:]
        fptr.close()
    
    if subset is None:
        subset = range(len(peak_list))

    bad_particles=None
    if bad_particle_list:
        if os.path.exists(bad_particle_list):
            bad_ptr = h5py.File(bad_particle_list,'r')
            bad_particles = bad_ptr['bad_particles'].value
            bad_ptr.close()
    leaf_indices={}

    min_dx = ds.index.get_smallest_dx()
    for clump in subset:
        this_clump = peak_list[clump]
        if not hasattr(this_clump,'units'):
            this_clump = ds.arr(this_clump,'code_length')

        region = ds.region(center     = this_clump,
                           left_edge  = this_clump-peak_radius*min_dx, 
                           right_edge = this_clump+peak_radius*min_dx)
        indices = region['particle_index'].astype('int64')
        if bad_particles is not None:
            mask = np.ones_like(indices,dtype='bool')
            for ni,i in enumerate(indices):
                #there's certainly a better way to do this than a loop.
                if i in bad_particles:
                    mask[ni]=False
            indices=indices[mask]
                    
        leaf_indices[clump]=indices
    return leaf_indices
    #for nc,indices in enumerate(leaf_indices):
     #   pw_full.annotate_select_particles(1.0, col='r', indices=indices)
   # pw_full.save(fname)
def shift_particles(ds=None, position=None,shift = np.zeros(3),shiftRight = False,grid_quan=None):
    """Shifts a periodically separated clump by the domain width.
    Looks for gaps in the positions larger than max('dx'), shifts one group
    to the right (left if shiftRight=Flase) to be spatially contiguous."""
    #max_dx may not be computed in the most efficient way.
    if ds is not None:  
        DomainLeft = ds.domain_left_edge
        DomainRight =  ds.domain_right_edge
        DomainWidth = DomainRight - DomainLeft
        max_dx = ds.index.grids[0].dds.max().in_units('code_length')
        min_dx = ds.index.get_smallest_dx()
    if grid_quan is not None:
        DomainLeft  = grid_quan["DomainLeft"]
        DomainRight = grid_quan["DomainRight"]
        DomainWidth = grid_quan["DomainWidth"]
        max_dx      = grid_quan["max_dx"]
        min_dx      = grid_quan["min_dx"]

    shifted=copy.copy(position)
    for i,axis in enumerate(['x','y','z']):

        dx = 'd'+axis
        nique = np.unique(shifted[:,i])
        nique.sort()

        #has to be close to the edges, or 'Periodic Wrap' isn't the problem.
        if np.abs(nique.max() - DomainRight[i]) > 3*max_dx:
            continue
        if np.abs(nique.min() - DomainLeft[i]) > 3*max_dx:
            continue
        delta_x = nique[1:] - nique[0:-1]
        max_delta_x = delta_x.max()
        break_index = np.where(delta_x == max_delta_x)
        if len(break_index) > 1:
            break_index = break_index[0]
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
            print( "indices inside", indices_late.sum())
        else:
            print( "indices inside", indices_late)
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


class NewSelectParticleCallback(PlotCallback):
    """
    Adds particle positions, based on a thick slab along *axis* with a
    *width* along the line of sight.  *p_size* controls the number of
    pixels per particle, and *col* governs the color.  *ptype* will
    restrict plotted particles to only those that are of a given type.
    *alpha* determines the opacity of the marker symbol used in the scatter.
    An alternate data source can be specified with *data_source*, but by
    default the plot's data source will be queried.
    """
    _type_name = "select_particles"
    region = None
    _descriptor = None
    _supported_geometries = ("cartesian", "spectral_cube", "cylindrical-2d")
    def __init__(self, width,indices=None, p_size=1.0, col='k', marker='o', stride=1,
                 ptype='all', minimum_mass=None, alpha=1.0, data_source=None):
        PlotCallback.__init__(self)
        self.width = width
        self.p_size = p_size
        self.color = col
        self.marker = marker
        self.stride = stride
        self.ptype = ptype
        self.minimum_mass = minimum_mass
        self.alpha = alpha
        self.data_source=data_source
        self.indices=indices
        if self.minimum_mass is not None:
            warnings.warn("The minimum_mass keyword is deprecated.  Please use "
                          "an appropriate particle filter and the ptype keyword instead.")


    def __call__(self, plot):
        data = plot.data
        if iterable(self.width):
            validate_width_tuple(self.width)
            self.width = plot.data.ds.quan(self.width[0], self.width[1])
        elif isinstance(self.width, YTQuantity):
            self.width = plot.data.ds.quan(self.width.value, self.width.units)
        else:
            self.width = plot.data.ds.quan(self.width, "code_length")
        # we construct a rectangular prism
        x0 = plot.xlim[0].to("code_length")
        x1 = plot.xlim[1].to("code_length")
        y0 = plot.ylim[0].to("code_length")
        y1 = plot.ylim[1].to("code_length")
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        if type(self.data_source)==YTCutRegion:
            mylog.warn("Parameter 'width' is ignored in annotate_particles if the "
                       "data_source is a cut_region. "
                       "See https://github.com/yt-project/yt/issues/1933 for further details.")
            self.region=self.data_source
        else:
            self.region=self._get_region((x0,x1), (y0,y1), plot.data.axis, data)
        ax = data.axis
        xax = plot.data.ds.coordinates.x_axis[ax]
        yax = plot.data.ds.coordinates.y_axis[ax]
        axis_names = plot.data.ds.coordinates.axis_name
        field_x = "particle_position_%s" % axis_names[xax]
        field_y = "particle_position_%s" % axis_names[yax]
        pt = self.ptype
        self.periodic_x = plot.data.ds.periodicity[xax]
        self.periodic_y = plot.data.ds.periodicity[yax]
        self.LE = plot.data.ds.domain_left_edge[xax], \
                  plot.data.ds.domain_left_edge[yax]
        self.RE = plot.data.ds.domain_right_edge[xax], \
                  plot.data.ds.domain_right_edge[yax]
        period_x = plot.data.ds.domain_width[xax]
        period_y = plot.data.ds.domain_width[yax]
        particle_x, particle_y = self._enforce_periodic(self.region[pt, field_x],
                                                        self.region[pt, field_y],
                                                        x0, x1, period_x,
                                                        y0, y1, period_y)
        gg = ( ( particle_x >= x0 ) & ( particle_x <= x1 )
           &   ( particle_y >= y0 ) & ( particle_y <= y1 ) )
        if self.indices is not None:
            mask_to_get = na.zeros(self.indices.shape, dtype='int32')
            found_any, mask = particle_ops.mask_particles(
                self.indices.astype('int64'), self.region['particle_index'].astype('int64'), mask_to_get)
            gg = ( gg & (mask == 1) )

        print( "nparticles in particle callback.", mask_to_get.sum())
        if self.minimum_mass is not None:
            gg &= (self.region[pt, "particle_mass"] >= self.minimum_mass)
            if gg.sum() == 0: return
        px, py = [particle_x[gg][::self.stride], particle_y[gg][::self.stride]]
        px, py = self._convert_to_plot(plot, [px, py])
        plot._axes.scatter(px, py, edgecolors='None', marker=self.marker,
                           s=self.p_size, c=self.color,alpha=self.alpha)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)

    def _enforce_periodic(self,
                          particle_x,
                          particle_y,
                          x0, x1, period_x,
                          y0, y1, period_y):
        #  duplicate particles if periodic in that direction AND if the plot
        #  extends outside the domain boundaries.
        if self.periodic_x and x0 > self.RE[0]:
            particle_x = uhstack((particle_x, particle_x + period_x))
            particle_y = uhstack((particle_y, particle_y))
        if self.periodic_x and x1 < self.LE[0]:
            particle_x = uhstack((particle_x, particle_x - period_x))
            particle_y = uhstack((particle_y, particle_y))
        if self.periodic_y and y0 > self.RE[1]:
            particle_y = uhstack((particle_y, particle_y + period_y))
            particle_x = uhstack((particle_x, particle_x))
        if self.periodic_y and y1 < self.LE[1]:
            particle_y = uhstack((particle_y, particle_y - period_y))
            particle_x = uhstack((particle_x, particle_x))
        return particle_x, particle_y

    def _get_region(self, xlim, ylim, axis, data):
        LE, RE = [None]*3, [None]*3
        ds = data.ds
        xax = ds.coordinates.x_axis[axis]
        yax = ds.coordinates.y_axis[axis]
        zax = axis
        LE[xax], RE[xax] = xlim
        LE[yax], RE[yax] = ylim
        LE[zax] = data.center[zax] - self.width*0.5
        LE[zax].convert_to_units("code_length")
        RE[zax] = LE[zax] + self.width
        if self.region is not None \
            and np.all(self.region.left_edge <= LE) \
            and np.all(self.region.right_edge >= RE):
            return self.region
        self.region = data.ds.region(data.center, LE, RE, data_source=self.data_source)
        return self.region

