from starter2 import *
#import three_loopers_1tff as tl

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

def return_core_count(self,n_min=1):
    output=[]
    for core_id in self.target_indices:
        if self.target_indices[core_id].size >= n_min:
            output.append(core_id)
    return output


class means_etc():
    def __init__(self,this_looper,n_min=3):
        thtr=this_looper.tr
        self.core_list = np.unique(this_looper.tr.core_ids)
        core_list=self.core_list
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
    fig=plt.figure(1, figsize=(8, 8))

    axScatter = plt.axes(rect_scatter)
    axHistx =   plt.axes(rect_histx)
    axHisty =   plt.axes(rect_histy)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    return fig, axScatter,axHistx, axHisty

