import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colorbar as cb
import matplotlib.colors as colors
import numpy as np
import yt
from importlib import reload

plt.clf()
plt.close('all')

if 1:
    prefix='ca03' #something to note the simulations you used.
    directory = '/scratch1/dcollins/Paper49_EBQU/ca03_turb_weak'
    frame = 339 
if 1:
    prefix='ca02' #something to note the simulations you used.
    directory = '/scratch1/dcollins/Paper49_EBQU/ca02_turb'
    frame=100
plot_directory = "%s/PlotsToTransfer"%os.environ['HOME']


#
# Create yt derived quantities.
#
def costhetaz(field,data):
    return data['magnetic_field_z']/data['magnetic_field_strength']
def thetaz(field,data):
    #arccos is behaving strangely with YTarrays
    #The zeros_like ensures we get the right type of object.
    output = np.zeros_like( data['costhetaz'])
    output[:] = np.arccos(data['costhetaz'][:].v) 
    return output

#This allows us to add these fields to any dataset

def add_cos_theta(obj):
    obj.add_field('costhetaz',costhetaz,units='dimensionless')
    obj.add_field('thetaz',thetaz,units='dimensionless')


#This conditional tests to see if region has been defined.
#if not, define it.  This process is slow so don't repeat if you don't have to.
ds = yt.load('%s/DD%04d/data%04d'%(directory,frame,frame))
add_cos_theta(ds)
region = ds.all_data()

#set up bins.
Nbins=64
theta_bins=np.linspace(0,np.pi,2*Nbins) #no particular reason other than debugging to use 2*Nbins
field_min = region['magnetic_field_strength'].min()
field_max = region['magnetic_field_strength'].max()
field_bins=np.logspace(np.log10(field_min),np.log10(field_max),Nbins)
bins={'thetaz':theta_bins,'magnetic_field_strength':field_bins,'magnetic_field_z':field_bins}

#
# Produce joint and marginalized distributions.
# Formally these are not PDFs since we don't explicityly normalize,
# BUT since the total volume is 1, they are in fact PDFs.
#
bin_fields = ['magnetic_field_strength','thetaz']
joint =     yt.create_profile(region, bin_fields=bin_fields,fields=['cell_volume'], weight_field=None, override_bins=bins)
prof_mag  = yt.create_profile(region, bin_fields=['magnetic_field_strength'],fields=['cell_volume'], weight_field=None, override_bins=bins)
prof_theta= yt.create_profile(region, bin_fields=['thetaz'],fields=['cell_volume'], weight_field=None, override_bins=bins)
prof_bz  =  yt.create_profile(region, bin_fields=['magnetic_field_z'],fields=['cell_volume'], weight_field=None, override_bins=bins)


#
# Now make plots with the profiles.
# 
#

if 1:
    #Marginalized distributions
    plt.close('all')
    fig,ax=plt.subplots(1,2)
    ax0=ax[0];ax1=ax[1]

    #Marginalize
    ph = joint['cell_volume']
    margin_theta= np.sum(ph,axis=0) #marginalize over magnetic field
    margin_mag = np.sum(ph,axis=1)

    thetabins = prof_theta.x_bins
    thetacen = 0.5*(thetabins[1:]+thetabins[:-1])

    #
    ax0.plot(thetacen,prof_theta['cell_volume'],c='k',label=r'$P(\theta)$')
    ax0.plot(thetacen,margin_theta,c='r',label= r'$\int P(B,\theta) dB$')
    ax0.set_title('theta')
    #ax0.legend(loc=2)

    b_bins = prof_bz.x_bins
    b_cen = 0.5*(b_bins[1:]+b_bins[:-1])

    ax1.plot(b_cen,prof_mag['cell_volume'],'k',label=r'$P(B)$')
    ax1.plot(b_cen,margin_mag,'r',label=r'$\int P(B,\theta) d\theta')
    ax1.set_title('B')
    ax1.set_xscale('log')

    plt.savefig('%s/%s_marginalized_distributions_%04d.png'%(plot_directory,prefix,frame))

        

#
# Check that P(B,theta) = P(B) P(theta)
#
if 1:
    #Joint and marginalized color plots
    fig2, axes2 = plt.subplots(1,2)
    ax0=axes2[0]; ax1=axes2[1] #this is in case I change the number of plots in the figure.
    p_theta=prof_theta['cell_volume']
    p_mag = prof_mag['cell_volume']

    #The product joint distribution P(B)P(theta)
    p_theta.shape = (p_theta.size,1)
    ph2 = p_theta*p_mag

    ph=np.transpose(joint['cell_volume'])
    min_val = ph[ph>0].min()
    max_val = ph.max()
    norm = colors.Normalize(vmin=min_val,vmax=max_val)
    extents = [prof_mag.x_bins.min(), prof_mag.x_bins.max(),prof_theta.x_bins.min(), prof_theta.x_bins.max()]
    pl=ax0.imshow(ph, norm=norm, interpolation='nearest',origin='lower',extent=extents,aspect='auto')
    cb=fig2.colorbar(pl,ax=ax0)
    cb.cmap.set_under('w')

    norm = colors.Normalize(vmin=min_val,vmax=max_val)
    pl1=ax1.imshow(ph2, norm=norm, interpolation='nearest',origin='lower', extent=extents,aspect='auto')
    cb1=fig2.colorbar(pl1,ax=ax1)
    cb1.cmap.set_under('w')

    fig2.savefig("%s/%s_joint_%04d.pdf"%(plot_directory,prefix,frame))
    plt.close(fig2)


if 1:
    #Joint and marginalized contour plots.
    #This works better than it has any right to.
    fig3, ax3 = plt.subplots(1,1)

    #I belive I may have a normalization wrong or something...
    # 
    ph2 = ph2/ph2.max()*ph.max()


    levels1 = np.arange(0,ph.max() , ph.max()/7)
    contour_sep=ax3.contour(ph2,levels1,colors='r',extent=extents)
    contour_joint=ax3.contour(ph ,levels1,colors='k',extent=extents)
    ax3.set_xlabel(r'$B$')
    ax3.set_ylabel(r'$\theta$')
    ax3.text(0.2,3.0,r'$P(B)P(\theta)$',color='r')
    ax3.text(0.2,2.7,r'$P(B,\theta)$',color='k')
    fig3.savefig("%s/%s_contours_%04d.pdf"%(plot_directory,prefix,frame))

#
# FINALLY we will predict P(B) from Pz(Bz)
# It doesn't work yet.
#

if 1:
    plt.close('all')
    fig,ax2=plt.subplots(1,1)

    dbin = b_cen[1:]-b_cen[:-1]
    pbz=prof_bz['cell_volume']
    dpbz=(pbz[1:]-pbz[:-1])/dbin
    ax2.plot(b_cen, b_cen*pbz)
    ax2.plot(b_cen[1:], -b_cen[1:]*dpbz,c='r',label= r'$-B dP/dz$')
    ax2.plot(b_cen, prof_mag['cell_volume'],c='k',label=r'$P(B)$')
    ax2.set_xscale('log')
    fig.savefig('%s/%s_prediction_%04d.png'%(plot_directory,prefix,frame))

