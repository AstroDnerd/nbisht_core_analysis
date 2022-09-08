from starter2 import *
import xtra_energy
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
sim_list=['u501']

from collections import defaultdict

class R_INFLECTION():
    def __init__(self, this_looper):
        self.this_looper=this_looper
        self.rinflection={}
        self.rinflection_list=[]

    def run(self,core_list=None, do_plots=False):
        this_looper=self.this_looper
        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        frame = this_looper.target_frame
        ds = this_looper.load(frame)
        G = ds['GravitationalConstant']/(4*np.pi)
        xtra_energy.add_energies(ds)
        for core_id in core_list:
            print('Potential %s %d'%(this_looper.sim_name,core_id))

            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            c = nar([ms.mean_x[-1], ms.mean_y[-1],ms.mean_z[-1]])
            
            R_SPHERE = 8/128
            rsph = ds.arr(R_SPHERE,'code_length')
            sp = ds.sphere(c,rsph)

            GE = np.abs(sp[YT_grav_energy_2])
            dv = np.abs(sp[YT_cell_volume])
            RR = sp[YT_radius]

            #2d distribution of GE vs r

            gbins = np.geomspace( GE[GE>0].min(), GE.max(),65)
            rbins = np.geomspace( RR [RR >0].min(), RR .max(),67)
            r_cen = 0.5*(rbins[1:]+rbins[:-1]) #we'll need this later.
            hist, xb, yb = np.histogram2d( RR , GE, bins=[rbins,gbins],weights=dv)


            #clear out annoying stragglers in the distribution.
            #any point that doesn't have any neighbors is eliminated.
            #h2 is the histogram to do math with.

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

                #the upper bound of the distribution
                #compute upper_envelope
                y = np.arange(h2.shape[1])
                y2d = np.stack([y]*h2.shape[0])
                argmax = np.argmax(y2d*(h2>0),axis=1)
                upper_envelope = gbins[argmax]
                keepers=upper_envelope>1

                #smooth for wiggles.
                UE = gaussian_filter(upper_envelope[keepers],1)

                #
                # Find the inflection point where the slope comes up again.
                #
                #the slope
                DUE = UE[1:]-UE[:-1]
                #the max up to this radius
                cummax=np.maximum.accumulate( UE)
                #where the slope is positive
                ok = DUE > 0
                #and not too close to the center (for noise)
                ok = np.logical_and(ok , r_cen[keepers][1:]>1e-3)
                #and it has to come down below half its maximum
                ok = np.logical_and(ok, UE[1:]<0.5*cummax[1:])

                index = np.argmax(r_cen[keepers])
                if ok.any():
                    #find out the radius where the inflection happens.
                    index = np.where(ok)[0][0]
                R_KEEP = r_cen[keepers][index]
                self.rinflection[core_id] = R_KEEP
                self.rinflection_list.append(R_KEEP)
            if do_plots:
                fig, ax = plt.subplots(1,1)
                pch.helper(hist,xb,yb,ax=ax,transpose=False)
                ax.plot(r_cen[keepers], UE, c='k')
                ax.scatter( r_cen[keepers][index], UE[index], c='r')
                axbonk(ax,xscale='log',yscale='log',xlabel='r',ylabel='rho')
                fig.savefig('plots_to_sort/r_inflection_%s_c%04d.png'%(this_looper.sim_name, core_id))

        return self.rinflection_list


if 0:
    import three_loopers_six as TL
    if 'inflection' not in dir():
        inflection = {}
        for sim in TL.loops:
            inflection[sim]=R_INFLECTION( TL.loops[sim])
        inflection[sim].run()
