from starter2 import *
import unyt

if 0:
    density_units=1000 #cm^-3
    mass_units_msun = 5900 #Msun
    length_units_au = 928000 #AU
    length_units_pc = 4.6 #pc
    velocity_units_km_s = 0.2
if 1:
    density_units=817 #cm^-3
    mass_units_msun = 2267 #Msun
    length_units_au = (3.6*unyt.pc).in_units('AU').v
    length_units_pc = 3.6 #pc
    velocity_units_km_s = 0.2
    time_units_Myr = 1.2
    time_units_s = (time_units_Myr*unyt.Myr).in_units('s').v

u_den = density_units*unyt.cm**-3
u_mass = mass_units_msun*unyt.Msun
u_length = length_units_au*unyt.AU
u_vel = velocity_units_km_s*unyt.km/unyt.s
u_eng = (u_mass*u_vel**2).in_units('erg')

color={'u05':'r','u10':'g','u11':'b'}
color.update({'u201':'r','u202':'g','u203':'b'})
color.update({'u301':'r','u302':'g','u303':'b'})
color.update({'u401':'r','u402':'g','u403':'b'})
color.update({'u501':'r','u502':'g','u503':'b'})
color.update({'u601':'r','u602':'g','u603':'b'})

mean_field={'u05':11.21, 'u10':3.545, 'u11':1.121}
mean_field.update({'u501':11.21, 'u502':3.545, 'u503':1.121})
mean_field.update({'u601':11.21, 'u602':3.545, 'u603':1.121})

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

def make_core_cmap(core_list, cmap='Spectral',seed=2653417):
    ncolors = max([ len(core_list), 400])
    if seed > 0:
        np.random.seed(seed)
        rands = np.random.random( len(core_list) ) * ncolors
        rands = np.round(rands)
    else:
        rands = np.linspace(0,ncolors, len(core_list)).astype('float')

    mymap = rainbow_map( ncolors, cmap = cmap)
    #mymap = rainbow_map( ncolors, cmap = 'tab20')
    color = [mymap(R) for R in rands]
    dic = dict( zip( core_list,color))
    return dic

G = 1620/(4*np.pi)
tff = np.sqrt( 3*np.pi/(32*G*1))
