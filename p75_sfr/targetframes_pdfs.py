
from starter2 import *
import track_loader as TL
import cfpack as cfp
from cfpack import stop,print
import testing.early_mask as em
import latex
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
reload(em)
if 0:
    np.set_printoptions(threshold=sys.maxsize)


class targetpdfs(): 
    def __init__(self,the_loop):
        print('init complete')
        self.this_looper = the_loop

    def toplot(prof,quan =YT_cell_volume): 
        xbins = prof.x_bins
        bin_center = 0.5*(xbins[1:]+xbins[:-1]) 
        pdf = prof[quan]
        return xbins, bin_center, pdf

    def vstime(sim):
        thtr = self.this_looper.tr
        last_frame = thtr.frames[-1]
        ds_last = self.this_looper.load(last_frame)  
        ad_last = ds_last.all_data()
        check = ad_last['density']>10000
        print(sim)
        has_true=any(check)
        print('some gas above rho=10^4 ',has_true)
        count_true = np.sum(check)  #plot this vs time
        print('how many ',count_true)
        return

    def coresvstime(self):
        print('inside cores count!!')
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        corecount = len(all_cores)
        return corecount

    def linearfit(x, a, b):
        y= a*x + b
        return y


    def getthepdfs(self, sim):
        print('inside cores count!!')
        thtr = self.this_looper.tr

        all_cores = np.unique(thtr.core_ids)
        core_list = all_cores
        frame = thtr.frames[0]

        ds = self.this_looper.load(frame)  
        em.add_tracer_density(ds)  
        ad = ds.all_data()
        deposit_tuple = ("deposit","target_particle_volume") 
 
        all_target_indices = self.this_looper.target_indices
        ad.set_field_parameter('target_indices',all_target_indices)
        ad.set_field_parameter('mask_to_get',np.zeros_like(all_target_indices,dtype='int32'))
        bins=None


        prof_alldata = yt.create_profile(ad,bin_fields=['density'],fields=['cell_volume'],weight_field=None, override_bins=bins)    
        prof_coresdata = yt.create_profile(ad,bin_fields=['density'],fields=[deposit_tuple],weight_field=None, override_bins=bins) 
      
        xbins1, bcen1, pdf1 = targetpdfs.toplot(prof_alldata)
        xbins2, bcen2, pdf2 = targetpdfs.toplot(prof_coresdata, deposit_tuple[1])
        ratio = pdf2/pdf1

        rets = cfp.fit(targetpdfs.linearfit, np.log10(bcen1[ratio>0]), np.log10(ratio[ratio>0]))
        pl = rets.popt[0] #then also plot power law value vs time
        fit_ratio = 10**targetpdfs.linearfit(np.log10(bcen1[ratio>0].v), *rets.popt)

        # 3 PLOTS PER FRAME (all_data, core_data)
        if 0: 
            outname = "plots_to_sort/sfpdfstest_%d_%s.png"%(frame,sim)
            fig,ax=plt.subplots(1,1)
            ax.plot(bcen1,pdf1,c='k',linewidth=1.0)
            ax.plot(bcen2,pdf2,c='k',linewidth=1.0, linestyle='dashed')
            ax.plot(bcen1[ratio>0],ratio[ratio>0],c='g',alpha=0.7) 
            ax.plot(bcen1[ratio>0],fit_ratio,c='r',alpha=0.5)

            ax.set(xlabel='cell_volume',ylabel=r'$PDF(\rho)$',xscale='log',yscale='log',ylim=(1e-7,1e0))
            fig.savefig(outname)
            print(outname)
            plt.close(fig)

        return pl



# YOU ARE HERE
sims=['m0230', 'm0231', 'm0232', 'm0233', 'm0234', 'm0235', 'm0236', 'm0237', 'm0238', 'm0239', 'm0240', 'm0241', 'm0242',\
      'm0243', 'm0244', 'm0245', 'm0246', 'm0247', 'm0250', 'm0260', 'm0270', 'm0280', 'm02100', 'm02110']  
#sims=['m02110']
if 'targets' not in dir():
    targets=True
# TO GET DATA AND STORE
if 1:
    TL.load_tracks(sims)
    if 'targets2' not in dir() or targets:
        pl = []
        cnum = []
        an_x = [30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,50,60,70,80,100,110]
        for sim in sims:
            targets = targetpdfs(TL.loops[sim])
            if 0: 
                numofcores = targets.coresvstime()
                cnum.append(numofcores)
            if 1:
                pl_val = targets.getthepdfs(sim)
                pl.append(pl_val)
            if 0:
                targets.vstime(sim)
        if 1:
            fig,ax=plt.subplots(1,1)
            if 0:  #number of cores vs time frames
                ax.scatter(an_x, cnum, c='g')
                ax.set(xlabel='target frame',ylabel='number of cores')
                outname='plots_to_sort/numcores_targettimes'
            if 1:  #slope vs time frames
                ax.scatter(an_x, pl, c='b')
                ax.set(xlabel='target frame',ylabel=r'slope, $P(*|\rho)$')
                outname='plots_to_sort/slope_targettimes'
            fig.savefig(outname)
            print(outname)
            plt.close(fig)


