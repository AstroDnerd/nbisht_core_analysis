from starter2 import *
import data_locations as dl
from collections import defaultdict

import davetools
reload(davetools)

plt.close('all')


class mass_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.angle_mean=defaultdict(list)
        self.angle_disp=defaultdict(list)
        self.cores_used=[]
    def run(self,core_list=None):
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
            self.ms = ms
            if ms.nparticles < 10:
                continue
            print('go ', core_id)
            self.cores_used.append(core_id)
            self.times = thtr.times


            for nf,frame in enumerate(thtr.frames):
                mask = ms.compute_unique_mask(core_id, dx=1./2048,frame=nf)
                density = thtr.c([core_id],'density')[mask,nf]
                cell_volume = thtr.c([core_id],'cell_volume')[mask,nf]
                bx = thtr.c([core_id],'magnetic_field_x')[mask,nf]
                by = thtr.c([core_id],'magnetic_field_y')[mask,nf]
                bz = thtr.c([core_id],'magnetic_field_z')[mask,nf]
                vx = thtr.c([core_id],'velocity_x')[mask,nf]
                vy = thtr.c([core_id],'velocity_y')[mask,nf]
                vz = thtr.c([core_id],'velocity_z')[mask,nf]

                bb =np.sqrt(bx*bx+by*by+bz*bz)
                vv =np.sqrt(vx*vx+vy*vy+vz*vz)
                BdotV = bx*vx+by*vy+bz*vz
                costheta=BdotV/(bb*vv)

                mean_cos_theta = (density*cell_volume*costheta).sum()/(density*cell_volume).sum()
                disp_cos_theta = (density*cell_volume*(costheta-mean_cos_theta)**2).sum()/(density*cell_volume).sum()
                self.angle_mean[core_id].append(mean_cos_theta)
                self.angle_disp[core_id].append(disp_cos_theta)

import three_loopers_1tff as tl
if 'clobber' not in dir():
    clobber=True
if 'mass_tool1' not in dir() or clobber:
    mass_tool1=mass_tool(tl.looper1)
    mass_tool1.run()


if 'mass_tool2' not in dir() or clobber:
    mass_tool2=mass_tool(tl.looper2)
    mass_tool2.run()
if 'mass_tool3' not in dir() or clobber:
    mass_tool3=mass_tool(tl.looper3)
    mass_tool3.run()

#
# Heat Map with Sample Tracks
#
fig,ax=plt.subplots(1,3, figsize=(12,4))
axes=ax.flatten()
if 1:
    for nt,tool in enumerate([mass_tool1,mass_tool2,mass_tool3]):


        #set up variables.
        G=1620./(4*np.pi)
        tff_global = np.sqrt(3*np.pi/(32*G*1))
        ntimes = tool.times.size
        ncores = len(tool.cores_used)
        angles = np.zeros([ntimes,ncores])
        these_times = tool.times/tff_global
        outname='plots_to_sort/avg_alignment_heatmap.pdf'

        #Make the angles into a 2d array we can plot
        for ncore,core_id in enumerate(tool.cores_used):
            this_ang = tool.angle_mean[core_id]
            angles[:,ncore]= np.arccos(this_ang)*180/np.pi

        #plot a few of the tracks
        nc = len(tool.cores_used)
        take_a_few = ((nc-1)*np.random.random(10)).astype('int')
        for ncore,core_id in enumerate(nar(tool.cores_used)[take_a_few]):
            this_ang=angles[:,ncore]
            axes[nt].plot(these_times, this_ang,c=[0.5]*4)

        #Make a 2d histogram.
        #There are probably better ways to do this.
        angle_bins_edge = np.linspace(0,180,32)
        angle_bins_edge = np.linspace(0,180,32)
        xbins = these_times
        ybins = 0.5*(angle_bins_edge[1:]+angle_bins_edge[:-1])
        nx = len(xbins) ; ny=len(ybins)
        TheX = np.r_[(ny)*[xbins]].transpose()
        TheY = np.r_[(nx)*[ybins]]

        hist = np.zeros( [xbins.size,ybins.size])
        for ntime, time in enumerate(these_times):
            thishist,bins = np.histogram(angles[ntime,:],bins=angle_bins_edge)
            hist[ntime,:]=thishist

        #Plot the map
        #First set up the color map, and force it to be white in areas where theres no points
        cmap = copy.copy(mpl.cm.get_cmap("viridis"))
        cmap.set_under('w')
        minmin = hist[hist>0].min()

        #set up the colorbar
        norm = mpl.colors.LogNorm(vmin=1,vmax=33)
        ploot=axes[nt].pcolormesh(TheX, TheY, hist, cmap=cmap,norm=norm,shading='nearest')
        axbonk(axes[nt],ylabel=None,xlabel=r'$t/t_{\rm{ff}}$',yscale='linear', ylim=[0,180])
    axes[0].set_ylabel(r'$<\theta>$')
    fig.colorbar(ploot)
    fig.savefig(outname)
    print(outname)

#
# Heat Map with Sample Tracks
#
fig,ax=plt.subplots(1,3, figsize=(12,4))
axes=ax.flatten()
if 1:
    for nt,tool in enumerate([mass_tool1,mass_tool2,mass_tool3]):


        #set up variables.
        G=1620./(4*np.pi)
        tff_global = np.sqrt(3*np.pi/(32*G*1))
        ntimes = tool.times.size
        ncores = len(tool.cores_used)
        angles = np.zeros([ntimes,ncores])
        these_times = tool.times/tff_global
        outname='plots_to_sort/avg_alignment_heatmap.pdf'

        #Make the angles into a 2d array we can plot
        for ncore,core_id in enumerate(tool.cores_used):
            this_ang = tool.angle_mean[core_id]
            angles[:,ncore]= np.arccos(this_ang)*180/np.pi
