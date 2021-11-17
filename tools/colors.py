from starter2 import *

color={'u05':'r','u10':'g','u11':'b'}
color.update({'u201':'r','u202':'g','u203':'b'})
color.update({'u301':'r','u302':'g','u303':'b'})
color.update({'u401':'r','u402':'g','u403':'b'})

#temporary definition, should be redudnant to davetools.
class rainbow_map():
    def __init__(self,n, cmap='jet'):
        norm = mpl.colors.Normalize()
        norm.autoscale(np.arange(n))
        #cmap = mpl.cm.jet
        self.color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)
    def __call__(self,val,n_fields=0):
        this_value = self.color_map.to_rgba(val)
        if n_fields > 0:
            this_value = [this_value]*n_fields
        return this_value

def make_core_cmap(core_list):
    np.random.seed(2653417)
    ncolors = max([ len(core_list), 400])
    mymap = rainbow_map( ncolors, cmap = 'Spectral')
    #mymap = rainbow_map( ncolors, cmap = 'tab20')
    rands = np.random.random( len(core_list) ) * ncolors
    color = [mymap(R) for R in rands]
    dic = dict( zip( core_list,color))
    return dic

