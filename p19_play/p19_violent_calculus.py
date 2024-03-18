
def calculus(obj, tsing):
        core_id=obj.core_id

        fig,axes=plt.subplots(1,2,figsize=(12,12))
        ax0=axes[0]; ax1=axes[1]
        #ax0=axes[0][1]; ax1=axes[0][0]
        #ax2=axes[1][0]; ax3=axes[1][1]
        #ax0=axes[1];ax1=axes[0];ax2=axes[2]
        #fig.subplots_adjust(wspace=0)
        rcen = 0.5*(obj.r_bins[1:]+obj.r_bins[:-1])

        #

        XXX,YYY = np.meshgrid( obj.times.flatten(),rcen)
        #norm = mpl.colors.LogNorm( dMdTf[dMdTf>0].min(), dMdTf.mean())
        ok = ~np.isnan(obj.meanvr)
        maxmax = np.abs(obj.meanvr[ok]).max()
        norm = mpl.colors.Normalize( -maxmax,maxmax)
        #pdb.set_trace()
        cmap=copy.copy(mpl.cm.get_cmap("seismic"))
        cmap.set_bad([0.9]*3)

        if 0:
            norm_velocity = mpl.colors.Normalize(vmin=-1000,vmax=1000)
            plot=ax0.pcolormesh( XXX,YYY, obj.divv,  norm=norm_velocity, shading='nearest', cmap=cmap)
            plot2=ax1.pcolormesh( XXX,YYY, obj.meanvr,  norm=norm_velocity, shading='nearest', cmap=cmap)
        if 1:
            norm_flux = mpl.colors.Normalize(vmin=-2,vmax=2)
            plot=ax0.pcolormesh(XXX,YYY,obj.dMdTf, norm=norm_flux,shading='nearest',cmap=cmap)
            plot2=ax1.pcolormesh(XXX,YYY,obj.dMdT3d, norm=norm_flux,shading='nearest',cmap=cmap)
        fig.colorbar(plot, label=r'divv',ax=ax0)
        fig.colorbar(plot2, label=r'vr',ax=ax1)
        ax0.set(yscale='log', xlabel='t/tff')
        ax1.set(yscale='log', xlabel='t/tff', ylabel='R [AU]')
        fig.savefig('plots_to_sort/calculus.png')
