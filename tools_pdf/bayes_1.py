from starter2 import *
# for new ds.all_data() profile plot - 
import core_dump
def rs(arr):
    return arr.reshape(arr.shape + (1,)) 
def rs2(arr):
    return arr.reshape(arr.shape + (1,) + (1,)) 

import trackage
reload(trackage)
if 'pdf_0000' not in dir():
    track = trackage.track_manager(None)

    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
    ds0 = yt.load(directory+"/DD0000/data0000")
    ad0 = ds0.all_data()
    extrema = {'density':[5e-3,100]}
    n_bins=128

    track.read('all_large_0000_0125.h5')
    den0 = rs2(track['density'][:,0]) #track['density'].reshape(track['density'].shape + (1,)) 
    den1 = rs2(track['density'][:,1]) #track['density'].reshape(track['density'].shape + (1,)) 
    cell_v0= rs2(track['cell_volume'][:,0] )
    cell_v1 = rs2(track['cell_volume'][:,1] )
    data = {'d1':den1,'density':den0, 'my_vol':cell_v0}#, 'kinetic_energy':ke} 
    data['KE'] =  rs2(track['kinetic_energy'][:,0])
    data['ME'] =  rs2(track['magnetic_energy'][:,0])

    bbox = np.array([[0.,1.]]*3) 
    ds_2 = yt.load_uniform_grid(data, den0.shape, length_unit="cm", bbox=bbox)
    ad_2 = ds_2.all_data() 
    pdf_0000 = yt.create_profile(ad0,'density',"cell_volume",
                                 weight_field=None,fractional=False,
                                 extrema=extrema, n_bins=n_bins)
    extrema['density']
    pdf_2 = yt.create_profile(ad_2,'density','my_vol',
                                 weight_field=None,fractional=False,
                                 extrema=extrema, n_bins=n_bins)
if 1:


    #extrema['ME']=
    import matplotlib.colors as colors
    pdf_2 = yt.create_profile(ad_2,['density','ME'],'my_vol',
                                 weight_field=None,fractional=False) #extrema=extrema, n_bins=n_bins
    plt.clf()
    norm = colors.LogNorm(vmin=1e-7,vmax=1e3)
    #norm = colors.Normalize( vmin=-7, vmax=3)#(vmin=1e-7,vmax=1e3)
    field=pdf_2['my_vol']
    myplot=plt.imshow(field, interpolation='nearest',origin='lower',norm=norm)
    #cb=plt.colorbar(myplot)
    #cb.cmap.set_under('w')
    plt.savefig('derp.png')


    #pdf_0 = yt.create_profile(ad0,['density','magnetic_energy'],'my_vol',
    #                             weight_field=None,fractional=False)
    #                             #extrema=extrema, n_bins=n_bins)
#

if 0:
    plt.clf()
    plt.plot(pdf_0000.x,pdf_0000['cell_volume'],label= r'$p(\rho)$')
    plt.plot(pdf_2.x,pdf_2['my_vol'],label= r'$p(\rho|*)$')
    plt.legend(loc=0)
    plt.xlabel(r'$\rho/\rho_0$')
    plt.ylabel(r'$PDF(\rho), PDF(\rho|*)$')
    plt.xscale('log');plt.yscale('log')
    plt.savefig('bayes_0.png')

if 0:
    import scatter_fit
    plt.clf()
    Pstar = 1
    Prho = pdf_0000['cell_volume']
    Prho_given = pdf_2['my_vol']
    Pstar_given = Prho_given*Pstar/Prho
    rho = pdf_0000.x
    ok = Pstar_given > 0
    plt.plot(rho[ok],Pstar_given[ok])
    scatter_fit.scatter_fit(plt,rho, Pstar_given, fit_range=[0.1,30], plot_points=False)
    plt.xlabel(r'$\rho/\rho_0$')
    plt.ylabel(r'$PDF(*|\rho)$')
    plt.xscale('log');plt.yscale('log')
    plt.savefig('bayes_1b.png')
if 0:
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(pdf_0000.x,pdf_0000['cell_volume'],label= r'$p(\rho)$', color='k')
    ax1.plot(pdf_2.x,pdf_2['my_vol'],label= r'$p(\rho|*)$',color='r')
    ax1.set_xscale('log');ax1.set_yscale('log')

    Pstar = 1
    Prho = pdf_0000['cell_volume']
    Prho_given = pdf_2['my_vol']
    Pstar_given = Prho_given*Pstar/Prho
    rho = pdf_0000.x
    ok = Pstar_given > 0
    ax2.plot(rho[ok],Pstar_given[ok], color='b')
    ax2.set_xscale('log');ax2.set_yscale('log')
    fig.savefig('bayes_1c.png')
