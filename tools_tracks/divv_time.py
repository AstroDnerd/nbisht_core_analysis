from starter2 import *
import data_locations as dl
from collections import defaultdict

import davetools
reload(davetools)

plt.close('all')

def dcolormesh(ax,xbins, ybins, hist):

    nx = len(xbins) ; ny=len(ybins)
    TheX = np.r_[(ny)*[xbins]].transpose()
    TheY = np.r_[(nx)*[ybins]]
    cmap = copy.copy(mpl.cm.get_cmap("viridis"))
    cmap.set_under('w')
    minmin = hist[hist>0].min()
    norm = mpl.colors.LogNorm(vmin=minmin,vmax=hist.max())
    ploot=ax.pcolormesh(TheX, TheY, hist, cmap=cmap,norm=norm,shading='nearest')
    return ploot

class div_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.unique_mass=defaultdict(list)
        self.divv=defaultdict(list)
        self.cores_used=[]
    def run(self,core_list=None, divmin=None, nybins=100):
        dx=1./2048
        nx = 1./dx
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            if ms.nparticles < 10:
                continue
            print('go ', core_id)
            self.cores_used.append(core_id)
            G=1620./(4*np.pi)
            tff_global = np.sqrt(3*np.pi/(32*G*1))

            self.times = thtr.times/tff_global


            divv = thtr.c([core_id],'velocity_divergence')
            cell_volume = thtr.c([core_id],'cell_volume')

            if divmin == None:
                #self.srt = np.sort(divv,axis=1)[:,0]
                self.srt = np.sort(np.abs(divv),axis=0)[:,0]
                dprob = 1/nybins
                iprob = int( dprob * self.srt.size)
                divmin = self.srt[iprob]
                


            maxmax = np.abs(divv).max()
            maxmax = 1000

            bin1 = np.logspace(np.log10(divmin), np.log10(maxmax),nybins+1)
            bin2 = -bin1[::-1]
            self.divvbins = np.concatenate([bin2,bin1])
            self.divvbins = np.linspace(-maxmax,maxmax,nybins)
            self.hists = [np.histogram( divv[:,nt], bins=self.divvbins)[0] for nt,time in enumerate(self.times)]
            self.hists = np.array(self.hists)

            fig,ax = plt.subplots(1,1)
            cen = 0.5*(self.divvbins[1:]+self.divvbins[:-1])
            ploot=dcolormesh(ax, self.times, cen, self.hists)
            axbonk(ax,xlabel=r'$t/t_{\rm{ff}}$',ylabel=r'$\nabla\cdot v$')
            #ax.set_yscale('symlog',linthresh=10*divmin)
            fig.colorbar(ploot)

            take_a_few = ((ms.nparticles-1)*np.random.random(10)).astype('int')
            for npart, part in enumerate(take_a_few):
                ax.plot(self.times, divv[npart,:], c=[0.5]*4)
            ax.plot(self.times,[0]*self.times.size,c='k')
            #ax.plot(self.times,divv.mean(axis=0),c='k')
            unique_div=[]
            unique_std=[]
            for nf, time in enumerate(self.times):
                mask2 = ms.compute_unique_mask(core_id, dx=1./2048,frame=nf)
                unique_div.append((divv[mask2,nf]*cell_volume[mask2,nf]).sum()/cell_volume[mask2,nf].sum())
                unique_std.append( np.sqrt(( ((divv[mask2,nf]-unique_div[-1])**2)*cell_volume[mask2,nf]).sum()/cell_volume[mask2,nf].sum()))

            ax.plot(self.times,divv.mean(axis=0),c='k')
            #ax.errorbar(self.times,divv.mean(axis=0),yerr=divv.std(axis=0),c='k')
            ax.errorbar(self.times,unique_div, yerr=unique_std,c='k')



            fig.savefig('plots_to_sort/divv_time_%s_c%04d.png'%(self.this_looper.out_prefix, core_id))
            #plt.close(fig)

            ax.clear()
            rm = rainbow_map(len(self.times))
            for nt, time in enumerate(self.times):
                the_x = np.sort(divv[:,nt])
                the_y = np.linspace(0,1,the_x.size)
                ax.plot(the_x,the_y,c=rm(nt))
            
            ax.grid(True)
            ax.set_xscale('symlog',linthresh=100)
            fig.savefig('plots_to_sort/divv_cdf_%s_c%04d.png'%(self.this_looper.out_prefix, core_id))

            plt.close(fig)


            #self.h = np.histogram()


import three_loopers_1tff as tl
if 'clobber' not in dir():
    clobber=True
if 'div_tool1' not in dir() or clobber:
    div_tool1=div_tool(tl.looper1)
    core_list = [31]
    #core_list = None
    div_tool1.run(core_list=core_list, divmin=None,nybins=100)

