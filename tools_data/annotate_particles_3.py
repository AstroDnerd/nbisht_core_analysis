
import yt
import matplotlib.pyplot as plt
from yt.visualization.plot_modifications import *
import pyximport; pyximport.install()
import particle_ops
import particle_grid_mask
from scipy.spatial import ConvexHull
import h5py
import time
import numpy as na
import os
import pdb
import copy
class TheseParticles(PlotCallback):
    """
    Adds particle positions, based on a thick slab along *axis* with a
    *width* along the line of sight.  *p_size* controls the number of
    pixels per particle, and *col* governs the color.  *ptype* will
    restrict plotted particles to only those that are of a given type.
    *alpha* determines the opacity of the marker symbol used in the scatter.
    An alternate data source can be specified with *data_source*, but by
    default the plot's data source will be queried.
    """
    _type_name = "these_particles2"
    region = None
    _descriptor = None
    _supported_geometries = ("cartesian", "spectral_cube", "cylindrical-2d")
    def __init__(self, width,positions, p_size=1.0, col='k', marker='o', stride=1,
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
        self.positions=positions
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
        px = self.positions[:,xax]
        py = self.positions[:,yax]
        if 0:
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
    _type_name = "select_particles4"
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
        print('AAA')
        
        if self.periodic_y and y0 < self.LE[1]:
            particle_y[ particle_y > y1 ] -= period_y
        if self.periodic_y and y1 > self.RE[1]:
            particle_y[ particle_y > y1 ] -= period_y
        if self.periodic_x and x0 < self.LE[0]:
            particle_x[ particle_x > x1 ] -= period_x
        if self.periodic_x and x1 > self.RE[0]:
            particle_x[ particle_x > x1 ] -= period_x
        if self.periodic_x and x0 > self.RE[0]:
            print('AAA')
            particle_x = uhstack((particle_x, particle_x + period_x))
            particle_y = uhstack((particle_y, particle_y))
        if self.periodic_x and x1 < self.LE[0]:
            print('AAA')
            particle_x = uhstack((particle_x, particle_x - period_x))
            particle_y = uhstack((particle_y, particle_y))
        if self.periodic_y and y0 > self.RE[1]:
            print('AAA')
            particle_y = uhstack((particle_y, particle_y + period_y))
            particle_x = uhstack((particle_x, particle_x))
        if self.periodic_y and y1 < self.LE[1]:
            print('AAA')
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

class ConvexHullCallback(PlotCallback):
    """
    Adds particle positions, based on a thick slab along *axis* with a
    *width* along the line of sight.  *p_size* controls the number of
    pixels per particle, and *col* governs the color.  *ptype* will
    restrict plotted particles to only those that are of a given type.
    *alpha* determines the opacity of the marker symbol used in the scatter.
    An alternate data source can be specified with *data_source*, but by
    default the plot's data source will be queried.
    """
    _type_name = "convex_hull_1"
    region = None
    _descriptor = None
    _supported_geometries = ("cartesian", "spectral_cube", "cylindrical-2d")
    def __init__(self, width,points=None, p_size=1.0, col='k', marker='o', stride=1,
                 ptype='all', alpha=1.0, data_source=None):
        PlotCallback.__init__(self)
        self.width = width
        self.p_size = p_size
        self.color = col
        self.marker = marker
        self.stride = stride
        self.ptype = ptype
        self.alpha = alpha
        self.data_source=data_source
        self.points=points

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
        pt = self.ptype
        self.periodic_x = plot.data.ds.periodicity[xax]
        self.periodic_y = plot.data.ds.periodicity[yax]
        self.LE = plot.data.ds.domain_left_edge[xax], \
                  plot.data.ds.domain_left_edge[yax]
        self.RE = plot.data.ds.domain_right_edge[xax], \
                  plot.data.ds.domain_right_edge[yax]
        period_x = plot.data.ds.domain_width[xax]
        period_y = plot.data.ds.domain_width[yax]
        raw_x = self.points[xax]
        raw_y = self.points[yax]
        particle_x, particle_y = self._enforce_periodic(raw_x,
                                                        raw_y,
                                                        x0, x1, period_x,
                                                        y0, y1, period_y)
        gg = ( ( particle_x >= x0 ) & ( particle_x <= x1 )
           &   ( particle_y >= y0 ) & ( particle_y <= y1 ) )
        px, py = [particle_x[gg][::self.stride], particle_y[gg][::self.stride]]
        px, py = self._convert_to_plot(plot, [px, py])

        points_2d = np.array(list(zip(px,py)))
        hull_2d = ConvexHull(points_2d)
        vert_list = np.concatenate([hull_2d.vertices, hull_2d.vertices[0:1]])
        plot._axes.plot(points_2d[vert_list,0], points_2d[vert_list,1], 'r')

        #plot._axes.scatter(px, py, edgecolors='None', marker=self.marker,
        #                   s=self.p_size, c=self.color,alpha=self.alpha)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)

    def _enforce_periodic(self,
                          particle_x,
                          particle_y,
                          x0, x1, period_x,
                          y0, y1, period_y):
        #  duplicate particles if periodic in that direction AND if the plot
        #  extends outside the domain boundaries.
        print('AAA')
        
        if self.periodic_y and y0 < self.LE[1]:
            particle_y[ particle_y > y1 ] -= period_y
        if self.periodic_y and y1 > self.RE[1]:
            particle_y[ particle_y > y1 ] -= period_y
        if self.periodic_x and x0 < self.LE[0]:
            particle_x[ particle_x > x1 ] -= period_x
        if self.periodic_x and x1 > self.RE[0]:
            particle_x[ particle_x > x1 ] -= period_x
        if self.periodic_x and x0 > self.RE[0]:
            print('AAA')
            particle_x = uhstack((particle_x, particle_x + period_x))
            particle_y = uhstack((particle_y, particle_y))
        if self.periodic_x and x1 < self.LE[0]:
            print('AAA')
            particle_x = uhstack((particle_x, particle_x - period_x))
            particle_y = uhstack((particle_y, particle_y))
        if self.periodic_y and y0 > self.RE[1]:
            print('AAA')
            particle_y = uhstack((particle_y, particle_y + period_y))
            particle_x = uhstack((particle_x, particle_x))
        if self.periodic_y and y1 < self.LE[1]:
            print('AAA')
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
