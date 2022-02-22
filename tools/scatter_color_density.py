
        if 1:
            bx_min=-0.5; bx_max=NumberOfOverlap.max()+1.5
            dx=1
            by_min=0; by_max=1.01; dy=0.01; Ny=int(1/dy)+1
            bins_x = np.arange(bx_min,bx_max,dx)
            bins_y = np.linspace(by_min,by_max,Ny)
            dy_actual = bins_y[1]-bins_y[0]
            h, xb,yb = np.histogram2d(NumberOfOverlap, Fraction,bins=(bins_x,bins_y))
            norm = colors.LogNorm(vmin=h[h>0].min(),vmax=h.max())
            cmap=copy.copy(mpl.cm.get_cmap("viridis"))
            #cmap.set_under('w')
            ix = np.floor(NumberOfOverlap).astype('int')
            iy = np.floor(Fraction/dy_actual).astype('int')
            vals = h[ix,iy]
            print(h[0,0])
            #print(vals)
            c=cmap(norm(vals))
            norms=c.sum(axis=1)
            #print( bins_x[ ix] - NumberOfOverlap)
            #pdb.set_trace()


            ax2[ns].scatter(NumberOfOverlap,Fraction,cmap=cmap,norm=norm,c=vals,s=1)
