
from starter2 import *
import xtra_energy
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import pcolormesh_helper as pch
import colors
import movie_frames 
G = colors.G
plt.close('all')
import scipy.stats

import radial_sphere 
reload(radial_sphere)

def mass_moment(prof,outname):
    ncol = 2
    nrow = 2
    fig,axes=plt.subplots(nrow,ncol, figsize=(12,12))
    if nrow==1:
        axes=[axes]
    if ncol==1:
        axes=[axes]

    radius = prof.r_cen*colors.length_units_au
    XXX,YYY = np.meshgrid(prof.times,radius)
    cmap=copy.copy(mpl.cm.get_cmap("seismic"))
    cmap.set_bad([0.9]*3)
    cmap.set_over('b')
    cmap.set_under('g')


    radial_sphere.dIdt(prof)
    field = 'M'
    if field in prof.fields:

        MMM = prof.fields[field]
        ok = ~np.isnan(MMM)
        radius_to_cut_off = 1000
        ok2 = ~np.isnan(MMM[:,-1])
        index = np.where( prof.r_cen[ok2]*colors.length_units_au < radius_to_cut_off)[0].max()
        mass_at_the_end = MMM[ok2,-1][index]
        fiducial=mass_at_the_end
        Max=MMM[ok].max()
        minner_chicken_dinner = fiducial**2/Max
        norm_mass = mpl.colors.LogNorm(vmin=minner_chicken_dinner,vmax=Max)

        nr = 0; nc = 0; ax = axes[nr][nc]
        plot=ax.pcolormesh(XXX,YYY,MMM,cmap=cmap, norm=norm_mass,shading='nearest')
        ax.set(yscale='log')
        fig.colorbar(plot,label='Mass',ax=ax)

    field = 'I'
    if field in prof.fields:

        III = prof.fields[field]
        ok = ~np.isnan(III)
        nr = 0; nc = 1; ax = axes[nr][nc]
        norm_mom = mpl.colors.LogNorm(vmin=III[ok].min(),vmax=III[ok].max())
        plot=ax.pcolormesh(XXX,YYY,III,cmap=cmap, norm=norm_mom,shading='nearest')
        fig.colorbar(plot,label='I',ax=ax)
        ax.set(yscale='log')
    field = 'dI'
    if field in prof.fields:
        nr = 1; nc = 0; ax = axes[nr][nc]
        dI = prof.fields[field]
        dXXX,dYYY = np.meshgrid(prof.tcen,radius)
        ok = ~np.isnan(dI)
        maxmax = np.abs(dI[ok]).max()
        #norm_dI = mpl.colors.SymLogNorm(vmin=-maxmax,vmax=maxmax,linthresh=10)
        norm_dI = mpl.colors.Normalize(vmin=-maxmax,vmax=maxmax)
        #norm_dI = mpl.colors.SymLogNorm(vmin=-100,vmax=100,linthresh=10)
        import dtools.davetools as dt
        #norm_dI = dt.norm_extrema(dI[ok].flatten())
        #norm_dI = mpl.colors.Normalize(vmin=-100,vmax=100)
        #norm_dI = mpl.colors.Normalize(vmin=dI[ok].min(),vmax=dI[ok].max())
        plot=ax.pcolormesh(dXXX,dYYY,dI,cmap=cmap,norm=norm_dI, shading='nearest')
        fig.colorbar(plot,label='dI',ax=ax)
        ax.set(yscale='log')
    if 0:
        nr=2; nc=0; ax=axes[nr][nc]
        dI = prof.fields[field]
        import dtools.math.equal_probability_binner as epb
        epb.equal_prob(dI.flatten(),32,ax=ax)
    field = 'ddI'
    if field in prof.fields:
        nr = 1; nc = 1; ax = axes[nr][nc]
        dI = prof.fields[field]
        dXXX,dYYY = np.meshgrid(prof.ddtcen,radius)
        ok = ~np.isnan(dI)
        #norm_ddI = mpl.colors.Normalize(vmin=dI[ok].min(),vmax=dI[ok].max())
        maxmax= np.abs(dI[ok]).max()
        norm_ddI = mpl.colors.Normalize(vmin=-maxmax,vmax=maxmax)

        print(maxmax)
        #norm_ddI = mpl.colors.SymLogNorm( vmin = -maxmax, vmax=maxmax,linthresh=0.1)
        #norm_ddI = mpl.colors.SymLogNorm( vmin = -0.01, vmax=0.01,linthresh=1e-5)
        #norm_ddI = mpl.colors.Normalize( vmin = -maxmax, vmax=maxmax)
        #norm_ddI = mpl.colors.LogNorm( vmin = np.abs(dI[ok][dI[ok]>0]).min(), vmax=np.abs(dI[ok]).max())
        plot=ax.pcolormesh(dXXX,dYYY,dI,cmap=cmap,norm=norm_ddI, shading='nearest')
        import dtools.math.equal_probability_binner as epb
        fig2,ax2=plt.subplots(1,1)
        epb.equal_prob(dI[ok].flatten(),32,ax=ax2)

        fig2.savefig('%s/derp'%plot_dir)
        



        fig.colorbar(plot,label='ddI',ax=ax)
        ax.set(yscale='log')


    af = axes.flatten()
    for ax in af:
        ax.set(xlim=af[0].get_xlim(), xlabel='t/tff',ylabel='R')
        tsing = prof.mon.get_tsing( prof.core_id)
        tsung = prof.mon.get_tsung( prof.core_id)
        ax.axvline(tsing, c=[0.5]*4)
        ax.axvline(tsung, c=[0.5]*4)
    fig.tight_layout()
    fig.savefig(outname)

def energy(prof,outname):
    ncol = 2
    nrow = 2
    fig,axes=plt.subplots(nrow,ncol, figsize=(12,12))
    if nrow==1:
        axes=[axes]
    if ncol==1:
        axes=[axes]

    radius = prof.r_cen*colors.length_units_au
    XXX,YYY = np.meshgrid(prof.times,radius)
    seismic=copy.copy(mpl.cm.get_cmap("seismic"))
    seismic.set_bad([0.9]*3)
    seismic.set_over('b')
    seismic.set_under('g')
    blue_r=copy.copy(mpl.cm.get_cmap("Blues_r"))
    blue_r.set_bad([0.9]*3)


    def ploot_color_at_final(MMM,ax):
        ok = ~np.isnan(MMM)
        radius_to_cut_off = 1000
        ok2 = ~np.isnan(MMM[:,-1])
        index = np.where( prof.r_cen[ok2]*colors.length_units_au < radius_to_cut_off)[0].max()
        mass_at_the_end = MMM[ok2,-1][index]
        fiducial=mass_at_the_end
        Max=MMM[ok].max()
        minner_chicken_dinner = fiducial**2/Max
        norm_mass = mpl.colors.LogNorm(vmin=minner_chicken_dinner,vmax=Max)
        plot=ax.pcolormesh(XXX,YYY,MMM,cmap=seismic, norm=norm_mass,shading='nearest')
        ax.set(yscale='log')
        fig.colorbar(plot,label=field,ax=ax)
    def ploot(MMM,ax,cmap=seismic,norm=None):
        ok = ~np.isnan(MMM)
        if norm is None:
            norm= mpl.colors.LogNorm( vmin=0.01,vmax=100)
        plot=ax.pcolormesh(XXX,YYY,MMM,cmap=cmap, norm=norm,shading='nearest')
        ax.set(yscale='log')
        fig.colorbar(plot,label=field,ax=ax)


    field = 'ET'
    if field in prof.fields:
        nr = 0; nc = 0; ax = axes[nr][nc]
        Grav = np.abs(prof.fields[field])
        ploot_color_at_final(Grav,ax)
    field = 'EG'
    if field in prof.fields:
        nr = 0; nc = 1; ax = axes[nr][nc]
        MMM = np.abs(prof.fields[field])/Grav
        ploot(MMM,ax)
    field = 'EB'
    if field in prof.fields:
        nr = 1; nc = 0; ax = axes[nr][nc]
        MMM = np.abs(prof.fields[field])/Grav
        ploot(MMM,ax)
    field = 'EK'
    if field in prof.fields:
        nr = 1; nc = 1; ax = axes[nr][nc]
        MMM = np.abs(prof.fields[field])/Grav
        ploot(MMM,ax)
        #ok = ~np.isnan(MMM)
        #norm = mpl.colors.LogNorm(vmin=MMM[ok].min(),vmax=MMM[ok].max())
        #ploot(MMM,ax,cmap=blue_r,norm=norm)

    af = axes.flatten()
    for ax in af:
        ax.set(xlim=af[0].get_xlim(), xlabel='t/tff',ylabel='R')
        tsing = prof.mon.get_tsing( prof.core_id)
        tsung = prof.mon.get_tsung( prof.core_id)
        ax.axvline(tsing, c=[0.5]*4)
        ax.axvline(tsung, c=[0.5]*4)
    fig.tight_layout()
    fig.savefig(outname)

def surface(surf,prof,outname):
    ncol = 4
    nrow = 2
    fig,axes=plt.subplots(nrow,ncol, figsize=(12,12))
    if nrow==1:
        axes=[axes]
    if ncol==1:
        axes=[axes]

    radius = prof.r_cen*colors.length_units_au
    XXX,YYY = np.meshgrid(prof.times,radius)
    seismic=copy.copy(mpl.cm.get_cmap("seismic"))
    seismic.set_bad([0.9]*3)
    seismic.set_over('b')
    seismic.set_under('g')
    blue_r=copy.copy(mpl.cm.get_cmap("Blues_r"))
    blue_r.set_bad([0.9]*3)


    def ploot(MMM,ax,cmap=seismic,norm=None):
        ok = ~np.isnan(MMM)
        maxmax = np.abs(MMM[ok]).max()
        minmax = MMM[ok][MMM[ok]>0].min()
        if norm is None:
            #norm= mpl.colors.Normalize(vmin=-maxmax, vmax=maxmax)
            norm= mpl.colors.LogNorm(vmin=minmax, vmax=maxmax)
        plot=ax.pcolormesh(XXX,YYY,MMM,cmap=cmap, norm=norm,shading='nearest')
        ax.set(yscale='log')
        fig.colorbar(plot,label=field,ax=ax)


    field = 'SK'
    if field in surf.field_list:
        nr = 0; nc = 0; ax = axes[nr][nc]
        arr = surf.fields[field]
        ploot(arr,ax)
        #import dtools.math.equal_probability_binner as epb
        #epb.equal_prob(surf.flatten(),32,ax=ax)
        kinetic = prof.fields['EK']
        ax=axes[1][0]
        norm = mpl.colors.LogNorm(vmin=0.01,vmax=100)
        ploot(arr/kinetic,ax=ax,norm=norm)

    af = axes.flatten()
#    for ax in af:
#        ax.set(xlim=af[0].get_xlim(), xlabel='t/tff',ylabel='R')
#        tsing = prof.mon.get_tsing( prof.core_id)
#        tsung = prof.mon.get_tsung( prof.core_id)
#        ax.axvline(tsing, c=[0.5]*4)
#        ax.axvline(tsung, c=[0.5]*4)
    fig.tight_layout()
    fig.savefig(outname)

def four_way(prof,outname):
    ncol = 4
    nrow = 4
    fig,axes=plt.subplots(nrow,ncol, figsize=(12,12))
    if nrow==1:
        axes=[axes]
    if ncol==1:
        axes=[axes]

    radius = prof.r_cen*colors.length_units_au
    XXX,YYY = np.meshgrid(prof.times,radius)
    seismic=copy.copy(mpl.cm.get_cmap("seismic"))
    seismic.set_bad([0.9]*3)
    seismic.set_over('b')
    seismic.set_under('g')

    def ploot_color_at_final(MMM,ax):
        norm_mass = mpl.colors.LogNorm(vmin=1e-6,vmax=1e-2)
        plot=ax.pcolormesh(XXX,YYY,MMM,cmap=seismic, norm=norm_mass,shading='nearest')
        return plot
    def ploot(MMM,ax,cmap=seismic,norm=None):
        ok = ~np.isnan(MMM)
        if norm is None:
            norm= mpl.colors.LogNorm( vmin=0.01,vmax=100)
        plot=ax.pcolormesh(XXX,YYY,MMM,cmap=cmap, norm=norm,shading='nearest')
        return plot


    field_list = ['EG','EK','ET','EB']
    for n1, f1 in enumerate(field_list):
        for n2, f2 in enumerate(field_list):
            ax=axes[n1,n2]
            arr1 = np.abs(prof.fields[f1])
            if f1 == f2:
                plot=ploot_color_at_final(arr1,ax)
                fig.colorbar(plot,label=f1,ax=ax)
            else:
                arr2 = np.abs(prof.fields[f2])
                MMM = arr2/arr1
                plot=ploot(MMM,ax)
                fname = '%s/%s'%(f2,f1)
                fig.colorbar(plot,label=fname,ax=ax)

    af = axes.flatten()
    for ax in af:
        ax.set(xlim=af[0].get_xlim(), xlabel='t/tff',ylabel='R',yscale='log')
        tsing = prof.mon.get_tsing( prof.core_id)
        tsung = prof.mon.get_tsung( prof.core_id)
        ax.axvline(tsing, c=[0.5]*4)
        ax.axvline(tsung, c=[0.5]*4)
    fig.tight_layout()
    fig.savefig(outname)

