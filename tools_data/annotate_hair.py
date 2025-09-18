
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

class HairCallback(PlotCallback):
    """
    Add streamlines to any plot, using the *field_x* and *field_y*
    from the associated data, skipping every *factor* datapoints like
    'quiver'. *density* is the index of the amount of the streamlines.
    *field_color* is a field to be used to colormap the streamlines.
    If *display_threshold* is supplied, any streamline segments where
    *field_color* is less than the threshold will be removed by having
    their line width set to 0.
    """

    _type_name = "hair"
    _supported_geometries = ("cartesian", "spectral_cube", "polar", "cylindrical")

    def __init__(
        self,
        field_x,
        field_y,
        factor=16,
        density=1,
        field_color=None,
        display_threshold=None,
        plot_args=None,
        this_looper=None,
        frame=None,
        core_id=None
    ):
        PlotCallback.__init__(self)
        def_plot_args = {}
        self.this_looper=this_looper
        self.frame=frame
        self.core_id=core_id
        self.field_x = field_x
        self.field_y = field_y
        self.field_color = field_color
        self.factor = factor
        self.dens = density
        self.display_threshold = display_threshold
        if plot_args is None:
            plot_args = def_plot_args
        self.plot_args = plot_args

    def __call__(self, plot):
        import pdb
        bounds = self._physical_bounds(plot)
        xx0, xx1, yy0, yy1 = self._plot_bounds(plot)

        # We are feeding this size into the pixelizer, where it will properly
        # set it in reverse order
        nx = plot.image._A.shape[1] // self.factor
        ny = plot.image._A.shape[0] // self.factor

        pixX = plot.data.ds.coordinates.pixelize(
            plot.data.axis, plot.data, self.field_x, bounds, (nx, ny)
        )
        pixY = plot.data.ds.coordinates.pixelize(
            plot.data.axis, plot.data, self.field_y, bounds, (nx, ny)
        )
        if self.field_color:
            field_colors = plot.data.ds.coordinates.pixelize(
                plot.data.axis, plot.data, self.field_color, bounds, (nx, ny)
            )

            if self.display_threshold:

                mask = field_colors > self.display_threshold
                lwdefault = matplotlib.rcParams["lines.linewidth"]

                if "linewidth" in self.plot_args:
                    linewidth = self.plot_args["linewidth"]
                else:
                    linewidth = lwdefault

                try:
                    linewidth *= mask
                    self.plot_args["linewidth"] = linewidth
                except ValueError as e:
                    err_msg = (
                        "Error applying display threshold: linewidth"
                        + "must have shape ({}, {}) or be scalar"
                    )
                    err_msg = err_msg.format(nx, ny)
                    raise ValueError(err_msg) from e

        else:
            field_colors = None

        X, Y = (
            np.linspace(xx0, xx1, nx, endpoint=True),
            np.linspace(yy0, yy1, ny, endpoint=True),
        )
        streamplot_args = {
            "x": X,
            "y": Y,
            "u": pixX,
            "v": pixY,
            "density": self.dens,
            "color": field_colors,
        }
        streamplot_args.update(self.plot_args)
        plot._axes.streamplot(**streamplot_args)
        plot._axes.set_xlim(xx0, xx1)
        plot._axes.set_ylim(yy0, yy1)

        #new stuff

        import trackage
        obnoxious_counter=0
        try:
            xax =plot.data.ds.coordinates.x_axis[plot.data.axis]
            yax =plot.data.ds.coordinates.y_axis[plot.data.axis]
            xaf = 'xyz'[xax]
            yaf = 'xyz'[yax]
            obnoxious_counter=1

            core_id = self.core_id
            this_looper=self.this_looper
            thtr=this_looper.tr
            ms = trackage.mini_scrubber(thtr,core_id)
            ms.particle_pos(core_id)
            frame_ind = np.where(this_looper.tr.frames ==self.frame)[0][0]
            obnoxious_counter=2
            all_x,all_y,all_z=ms.particle_x,ms.particle_y, ms.particle_z
            all_p = [all_x,all_y,all_z]
            all_p_s = np.stack(all_p)
            max_max = all_p_s.max(axis=1).max(axis=1)
            min_min = all_p_s.min(axis=1).min(axis=1)
            cen = 0.5*(min_min+max_max)
            XX,YY= all_p[xax].transpose(), all_p[yax].transpose()
            plot._axes.scatter(XX[0,:],YY[0,:], c='k')
            plot._axes.plot(XX[:frame_ind+1,:],YY[:frame_ind+1,:], c=[0.5]*4, zorder=7, linewidth=0.1)
            plot._axes.scatter(XX[frame_ind,:],YY[frame_ind,:], c='r', s=1,zorder=1)
        except:
            print("Failed by ",obnoxious_counter)
            pdb.set_trace()
        #plot._axes.set_xlim(xx0, xx1)
        #plot._axes.set_ylim(yy0, yy1)




