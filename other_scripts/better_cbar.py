#
# Colorbars are a pain in the ass in matplotlib.  
# Too many fine grained objects one has to make.
# I want to make one object.
# 
import matplotlib
import numpy as na
import pdb

try:
    algae = matplotlib.cm.get_cmap("algae")
except:
    pass
class better_cbar(object):
    def __init__(self,values=None,zlim=None,cb_axes=None,scale='linear',cmap='jet',label='None',
                clear_axes_ticks=True, under_white=False):
        self.values=values
        self.zlim=zlim
        self.cb_axes=cb_axes
        self.scale=scale
        print 'better_cbar:', cmap
        self.cmap=cmap
        self.label=label
        if self.values == None and self.zlim == None:
            print "One must supply either zlim or value range"

        if self.scale=='linear':
            norm_object = matplotlib.colors.Normalize
        else:
            norm_object = matplotlib.colors.LogNorm
        if self.zlim:
            self.norm = norm_object(self.zlim[0],self.zlim[1])
        else:
            self.norm = norm_object()
            self.norm.autoscale(self.values)
        
        if self.cmap == 'algae':
            self.cmap = algae
        else:
            self.cmap = matplotlib.cm.__dict__[self.cmap]
        if under_white:
            self.cmap.set_under('w')
        self.color_map = matplotlib.cm.ScalarMappable(norm=self.norm,cmap=self.cmap)
        self.to_rgba = self.color_map.to_rgba
        if self.cb_axes is not None:
            self.cb_base=matplotlib.colorbar.ColorbarBase(cmap=self.cmap,norm=self.norm,ax=self.cb_axes)
            self.cb_axes.yaxis.set_label_position('right')
            self.cb_axes.set_ylabel(self.label)
        self.minval = None
        self.maxval = None
        self.values = []
    def map(self,values):
        self.values.append(values)
        if self.minval:
            self.minval = min(self.minval,min(values))
            self.maxval = max(self.maxval,max(values))
        else:
            try:
                self.minval = min(values)
                self.maxval = max(values)
            except:
                self.minval = values
                self.maxval = values
        try:
            len(values)
            return [self.to_rgba(v) for v in values] 
        except:
            return self.to_rgba(values)
    def extents(self):
        return self.minval,self.maxval
