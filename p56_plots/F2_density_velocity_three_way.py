
from starter2 import *
import matplotlib.image as mpimg

import data_locations as dl
reload(dl)
reload(trackage)
plt.close('all')

form='pdf'
import three_loopers_1tff as tl

tm = rainbow_map(15)
#all_cores = all_cores[:10]
def make_prof(ds,fields,weight_field=None,accumulation=False,fractional=True,n_bins=64,extrema=None):
    reg = ds.all_data()
    prof = yt.create_profile(reg,fields[0],fields[1] ,weight_field=weight_field,accumulation=accumulation,
                            fractional=fractional, n_bins=n_bins, extrema=extrema)
    the_x = 0.5*(prof.x_bins[1:]+prof.x_bins[0:-1])
    the_y = prof[fields[1]]
    #if units[0] is not None:
    #    the_x = the_x.in_units(units[0])
    #if units[1] is not None and fractional is not True:
    #    the_y = the_y.in_units(units[1])
    output={}
    output['prof']=prof
    output['the_x']=the_x
    output['the_y']=the_y
    return output

plt.clf()
odir=os.environ['HOME']+'/PigPen/'
odir = "./plots_to_sort/"


def return_core_count(self,n_min=1):
    output=[]
    for core_id in self.target_indices:
        if self.target_indices[core_id].size >= n_min:
            output.append(core_id)
    return output


class means_etc():
    def __init__(self,this_looper,n_min=3):
        thtr=this_looper.tr
        core_list = return_core_count(this_looper,n_min)
        self.dmeans = np.zeros_like(core_list,dtype='float')
        self.dstds = np.zeros_like(core_list,dtype='float')
        self.d_logmeans = np.zeros_like(core_list,dtype='float')
        self.d_logstds  = np.zeros_like(core_list,dtype='float')
        self.v_logmeans = np.zeros_like(core_list,dtype='float')
        self.v_logstds  = np.zeros_like(core_list,dtype='float')
        self.variance = np.zeros_like(core_list,dtype='float')
        self.vmeans    = np.zeros_like(core_list,dtype='float')
        self.vstds = np.zeros_like(core_list,dtype='float')
        self.npart = np.zeros_like(core_list,dtype='float')
        self.vrel  =  np.zeros_like(core_list,dtype='float')
        self.volume =  np.zeros_like(core_list,dtype='float')
        self.temp_all_rho=[]

        for i,nc in enumerate(core_list):
            this_density = thtr.c(int(nc),'density')[:,0]
            self.temp_all_rho.append(this_density)
            this_vel = thtr.c(int(nc),'velocity_magnitude')[:,0]
            this_volume = thtr.c(int(nc),'cell_volume')[:,0]
            self.npart[i] = this_density.size
            self.volume[i] = this_volume.sum()
            self.dmeans[i]=this_density.mean()
        
            self.dstds[i] = this_density.std()
            self.d_logmeans[i]=np.exp(np.log(this_density).mean())
            self.d_logstds[i] =np.exp(np.log(this_density).std())
            self.v_logmeans[i]=np.exp(np.log(this_vel).mean())
            self.v_logstds[i] =np.exp(np.log(this_vel).std())
            self.vmeans[i]= this_vel.mean()
            self.vstds[i] = this_vel.std()
            ms = trackage.mini_scrubber(thtr,[nc],do_velocity=True)
            self.variance[i] =  np.mean(ms.rel_vx**2)+np.mean(ms.rel_vy**2)+np.mean(ms.rel_vz**2)

    def make_profiles(self,looper):
        #note to future self, I put this here but didn't finish.
        if 'prof0' not in dir() and False:
            frame = dl.target_frames[this_simname]
            ds0 = yt.load("%s/DD%04d/data%04d"%(dl.sims[this_simname],0,0))
            ds1 = yt.load("%s/DD%04d/data%04d"%(dl.sims[this_simname],frame,frame))
            prof0=make_prof(ds0,['density','cell_volume'])
            prof1=make_prof(ds1,['density','cell_volume'])

        if 'prof_vel' not in dir():
            prof_vel=make_prof(ds1,['velocity_magnitude','cell_volume'])



m1 = means_etc(tl.looper1 )
m2 = means_etc(tl.looper2 )
m3 = means_etc(tl.looper3 )

tl.looper1.c='r'
tl.looper2.c='g'
tl.looper3.c='b'

def three_way_bean():
    from matplotlib.ticker import NullFormatter
    nullfmt = NullFormatter()         # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(1, figsize=(8, 8))

    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    return axScatter,axHistx, axHisty

#    # the scatter plot
#    axScatter.scatter(x, y)
#
#    # now determine nice limits by hand
#    binwidth = 0.25
#    xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
#    lim = (int(xymax/binwidth) + 1) * binwidth
#
#    axScatter.set_xlim((-lim, lim))
#    axScatter.set_ylim((-lim, lim))
#
#    bins = np.arange(-lim, lim + binwidth, binwidth)
#    axHistx.hist(x, bins=bins)
#    axHisty.hist(y, bins=bins, orientation='horizontal')
#
#    axHistx.set_xlim(axScatter.get_xlim())
#    axHisty.set_ylim(axScatter.get_ylim())
##
#    plt.savefig(outname)

   


if 1:
    if 1:
        ax, ax_den_hist,ax_vel_hist=three_way_bean()
        ok = slice(None)
        ax.scatter(np.log10(m1.dmeans[ok]),np.log10(m1.variance[ok])/2,c=tl.looper1.c,label=tl.looper1.out_prefix, s=0.1)
        ax.scatter(np.log10(m2.dmeans[ok]),np.log10(m2.variance[ok])/2,c=tl.looper2.c,label=tl.looper2.out_prefix, s=0.1)
        ax.scatter(np.log10(m3.dmeans[ok]),np.log10(m3.variance[ok])/2,c=tl.looper3.c,label=tl.looper3.out_prefix, s=0.1)


        ax_den_hist.hist(np.log10(m1.dmeans[ok]), histtype='step',color=tl.looper1.c,label=tl.looper1.out_prefix)
        ax_den_hist.hist(np.log10(m2.dmeans[ok]), histtype='step',color=tl.looper2.c,label=tl.looper2.out_prefix)
        ax_den_hist.hist(np.log10(m3.dmeans[ok]), histtype='step',color=tl.looper3.c,label=tl.looper3.out_prefix)

        ax_vel_hist.hist(np.log10(m1.variance[ok])/2, histtype='step', orientation='horizontal',color=tl.looper1.c)
        ax_vel_hist.hist(np.log10(m2.variance[ok])/2, histtype='step', orientation='horizontal',color=tl.looper2.c)
        ax_vel_hist.hist(np.log10(m3.variance[ok])/2, histtype='step', orientation='horizontal',color=tl.looper3.c)
        axbonk(ax,yscale='linear', xscale='linear',  ylabel=r'$\log_{10} \sigma_v$', xlabel=r'$\log_{10} \langle\rho\rangle$')
        axbonk(ax_vel_hist,yscale='linear', xscale='linear',  ylabel=None, xlabel=r'$N$')
        axbonk(ax_den_hist,yscale='linear', xscale='linear',  ylabel=r'$N$', xlabel=None)
        ax_den_hist.legend(loc=1)
        ax_den_hist.set_xticks([])
        ax_vel_hist.set_yticks([])
        ax.set_xlim(-1,2)
        ax_den_hist.set_xlim(ax.get_xlim())
        ax_vel_hist.set_ylim(ax.get_ylim())

        plt.savefig(odir+'/pre_rho_rms_v_rms.%s'%form)

