
from starter1 import *
import matplotlib.colors as colors
import yt
import looper
reload(looper)
from scipy.stats import kurtosis
from scipy.stats import skew
from matplotlib.ticker import LogLocator, LogFormatterSciNotation as LogFormatter
import trackage
reload(trackage)
def rs(arr):
    return arr.reshape(arr.shape + (1,)) 
def rs2(arr):
    return arr.reshape(arr.shape + (1,) + (1,)) 
if 'always' not in dir():
    always = False
def expform(float_in, format = "%0.1e"):
    """exponent format: 1.5e+03 goes to $1.5 \times 10^{3}$"""
    str1 = format%float_in
    tex_string = "$"
    try:
        exp_pos = str1.index("e")
        exponent_shift = 1
        if str1[exp_pos+1] == "+":
            exponent_shift+1
        mantissa = str1[0:exp_pos]
        exponent = str1[exp_pos+exponent_shift:]
        if float(mantissa) != 1:
            tex_string += mantissa + "\\times "
        tex_string += "10^{" 
        tex_string += str(int(exponent))
        tex_string += "}$"
    except:
        tex_string += str1
        tex_string += "$"
    return '%s'%tex_string
plt.close('all')
def tick_fixer(xticks,xbins):
    x_slice=slice(None)
    if xticks[-1] > len(xbins): x_slice=slice(None,-1)
    new_ticklabels = [ expform(xbins[n]) for n in xticks[x_slice]]
    return new_ticklabels

if 'track' not in dir() or always:
    track = trackage.track_manager(None)

    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
    #track.read('tracks_first_last_all_nonzero.h5')
    track.read('tracks_frame012_all_nonzero.h5')

    #fname_initial = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/DD0000/data0000'
    fname_initial = '/home/dcollins/scratch/u05-r4-l4-128/DD0000/data0000'
    fname_initial = '/home/dcollins/scratch/u05-r4-l4-128/DD0001/data0001'
    fname_initial = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/DD0001/data0001'

    ds_all = yt.load(fname_initial)
    ad_all = ds_all.all_data()

if 'ad_pre' not in dir() or always:
    mag = rs2(track['magnetic_field_strength'][:,1])
    vel = rs2(track['velocity_magnitude'][:,1])
    #vel_bulk = track['_vel_mag'].reshape(track['_vel_mag'].shape + (1,))
    #vel_div = track['velocity_divergence'].reshape(track['velocity_divergence'].shape + (1,))
    #div_vel = track['divvel_plus'].reshape(track['divvel_plus'].shape + (1,))
    den = rs2(track['density'][:,1])
    cell_v = rs2(track['cell_volume'][:,1])
     
    data = {'density':den, 
            'magnetic_field_strength':mag, 
            #'_vel_mag':vel_bulk, 
            #'velocity_divergence':vel_div, 
            'my_vol':cell_v, 
            #'divvel_plus':div_vel,
            'velocity_magnitude':vel}  
    bbox = np.array([[0.,1.]]*3)
    ds = yt.load_uniform_grid(data, den.shape, length_unit="cm", bbox=bbox)
    ad = ds.all_data() 

    extrema   = {'density':[5e-3,100],
                 'magnetic_field_strength':[0.25,100],
                 'velocity_magnitude':[5e-3,100]} 
def ratio_pair( field,  ad_pre=None, ad_cores=None,extrema={}, zlim=[1e-2,1e4]):
    pdf_pre = yt.create_profile(ad_pre, [field[0], field[1]], field[2], 
                                weight_field=None, fractional=False,
                                extrema=extrema)
    fig, ax = plt.subplots() 
    
    norm = colors.LogNorm(vmin=zlim[0],vmax=zlim[1])  # vmin=1e-4, vmax=1 
    norm = None
    z_pre = pdf_pre[field[2]]
    pre_plot = ax.pcolormesh(z_pre, norm=norm)  ##changed plt to ax   

    #pre_plot.set_norm(norm)
    #locator = LogLocator()
    #formatter = LogFormatter()
    #cbar = fig.colorbar(pre_plot,ax=ax,norm=norm)  ## can add ticks=[]
    #Set color to be used for low out-of-range values:
    #cbar.cmap.set_under('w')

    #cbar.locator = locator
    #cbar.formatter = formatter
    #cbar.update_normal(pre_plot) 
    #cbar.set_cmap('arbre')
    #cbar.set_clim(vmin=1e-2,vmax=1e4) 
    #cbar.ax.set_ylabel(r'density ($g/cm^3$)')  #(r'cell_volume ($cm^3$)') 

    ax.set_xlabel(r'%s'%(field[0]))
    ax.set_ylabel(r'%s'%(field[1])) 

    xticks=ax.get_xticks().astype('int')
    these_xbins = pdf_pre.x_bins
    new_labels = tick_fixer( xticks, these_xbins)
    ax.set_xticklabels(new_labels)

    yticks=ax.get_yticks().astype('int')
    these_ybins = pdf_pre.y_bins
    new_labels = tick_fixer( yticks, these_ybins)
    ax.set_yticklabels(new_labels,rotation='45')

    outname = '/home/dcollins4096/PigPen/plot_test_%s_%s_%s.png'%tuple(field)
    plt.savefig(outname)
    print(outname)
    return {'pdf_pre':pdf_pre,'ax':ax}

this_extrema={}
for field in ['density','velocity_magnitude']:
    this_extrema[field]=extrema[field]
stuff=ratio_pair(['density','velocity_magnitude','cell_volume'], 
           ad_pre=ad_all, ad_cores=ad,extrema=this_extrema)
