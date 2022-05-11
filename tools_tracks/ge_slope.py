from starter2 import *
import xtra_energy
import three_loopers_u500 as TL
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
sim_list=['u501']

def plot_phi(this_looper,core_list=None):
    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    frame = this_looper.target_frame
    for core_id in core_list:
        print('Potential %s %d'%(this_looper.sim_name,core_id))

        ds = this_looper.load(frame)
        xtra_energy.add_energies(ds)
        ms = trackage.mini_scrubber(this_looper.tr,core_id)
        c = nar([ms.mean_x[-1], ms.mean_y[-1],ms.mean_z[-1]])
        
        rsph = ds.arr(8.0/128,'code_length')
        sp = ds.sphere(c,rsph)

        GE = np.abs(sp['grav_energy'])
        dv = np.abs(sp['cell_volume'])
        RR = sp['radius']
        r_sphere=sp['radius']
        DD = sp['density']

        fig,ax=plt.subplots(1,1)
        ax0=ax

        #2d distribution of GE and r

        rbins = np.geomspace( r_sphere[r_sphere>0].min(), r_sphere.max(),67)
        r_cen = 0.5*(rbins[1:]+rbins[:-1]) #we'll need this later.
        gbins = np.geomspace( GE[GE>0].min(), GE.max(),65)
        hist, xb, yb = np.histogram2d( r_sphere, GE, bins=[rbins,gbins],weights=dv)

        #clear out annoying stragglers in the distribution.
        #any point that doesn't have any neighbors is eliminated.
        if 1:
            #h2 is the histogram.
            #we'll remove any stragglers.
            h2 = hist+0
            shifter = np.zeros(nar(h2.shape)+2)
            cuml = np.zeros(h2.shape)
            c_center = slice(1,-1)
            #hb is just "is h2 nonzero"
            #we'll slide this around to look for neighbors
            hb = (h2>0)
            shifter[c_center,c_center]=hb
            nx,ny=shifter.shape
            for i in [0,1,2]:
                for j in [0,1,2]:
                    if i==1 and j==1:
                        continue
                    s1 = slice(i,nx-2+i)
                    s2 = slice(j,ny-2+j)
                    cuml += shifter[s1,s2]
            #kill points that don't have neighbors.
            h2 *= (cuml >0)

        if 1:
            #Compute the upper bound of the histogram
            #smooth it
            #look for the point where the slope goes up.
            #but to avoid wiggles, it has to first come down a lot.
            y = np.arange(h2.shape[1])
            y2d = np.stack([y]*h2.shape[0])
            argmax = np.argmax(y2d*(h2>0),axis=1)
            upper_envelope = gbins[argmax]
            keepers=upper_envelope>1
            #smooth for wiggles.
            UE = gaussian_filter(upper_envelope[keepers],1)

            #the slope
            DUE = UE[1:]-UE[:-1]
            #the max up to this radius
            cummax=np.maximum.accumulate( UE)
            #where the slope is positive
            ok = DUE > 0
            #and not too close
            ok = np.logical_and(ok , r_cen[keepers][1:]>1e-3)
            #and it has to come down below half its maximum
            ok = np.logical_and(ok, UE[1:]<0.5*cummax[1:])


        if ok.any():
            #find out the radius where the inflection happens.
            index = np.where(ok)[0][0]
            R_KEEP = r_cen[keepers][index]

        if 1:
            #do the fit.
            ok_fit = np.logical_and(r_sphere < R_KEEP, GE>0)
            ok_fit = np.logical_and(ok, r_sphere > r_cen[keepers].min())
            rok=r_sphere[ok_fit]
            def powerlaw( x, q, r0):
                return q*x+r0
            popt, pcov=curve_fit(powerlaw, np.log10(rok), np.log10(GE[ok_fit]))
            #Dan, popt will have the fit parameters you want.

        if 1:
            #plot stuff
            ax.plot(r_cen[keepers], UE)
            ax.scatter( r_cen[keepers][index],UE[index],c='r')
            pch.helper(h2,xb,yb,ax=ax0,transpose=False)
            ax.plot( rok, 10**powerlaw(np.log10(rok), *popt),c='r')


        if 0:
            #color the upper envelope
            #to make sure we get it right.
            print(hist.shape)
            y = np.arange(hist.shape[1])
            y2d = np.stack([y]*hist.shape[0])
            argmax = np.argmax(y2d*(hist>0),axis=1)
            ind = np.arange( hist.shape[0])
            #take = np.ravel_multi_index(nar([argmax,ind]),hist.shape)
            take = np.ravel_multi_index(nar([ind,argmax]),hist.shape)
            h1=hist.flatten()
            h1[take]=hist.max()
            h1.shape=hist.shape
            pch.helper(h1,xb,yb,ax=ax0,transpose=False)
        outname='plots_to_sort/%s_c%04d_potfit'%(this_looper.sim_name,core_id)
        axbonk(ax0,xscale='log',yscale='log',xlabel='r',ylabel='grad phi sq')
        fig.savefig(outname)
        print(outname)

for sim in TL.loops:
    if sim == 'u501':
        continue
    all_cores=np.unique( TL.loops[sim].tr.core_ids)
    core_list=list(all_cores)
    plot_phi( TL.loops[sim])#,core_list=core_list)
