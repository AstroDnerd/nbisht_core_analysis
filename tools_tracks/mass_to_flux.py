from starter2 import *
import data_locations as dl
from collections import defaultdict

import davetools
reload(davetools)

plt.close('all')

min_r = 1./2048
G = 1620/(4*np.pi)

class mphi_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.Energy_G_core = defaultdict(list)
        self.Energy_B_core = defaultdict(list)
        self.Energy_M_core = defaultdict(list)
        self.Energy_K_core = defaultdict(list)
        self.Energy_W_core = defaultdict(list)
        self.cores_used  = []

        self.Energy_G_time = defaultdict(list)
        self.Energy_B_time = defaultdict(list)
        self.Energy_M_time = defaultdict(list)
        self.costheta      = defaultdict(list)
    def run(self,core_list=None, hull=None):
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

            Rhull = None
            if hull is not None:
                if core_id not in hull.cores_used:
                    continue
                else:
                    this_hull = hull.hulls[core_id]
                    Rhull = this_hull.volume**(1./3)
            self.cores_used.append(core_id)
            self.times = thtr.times


            all_density = thtr.c([core_id],'density')
            all_phi = thtr.c([core_id],'PotentialField')
            all_cell_volume = thtr.c([core_id],'cell_volume')
            all_magnetic_field_strength = thtr.c([core_id],'magnetic_energy')
            all_kinetic_energy = thtr.c([core_id],'kinetic_energy')
            all_vorticity = thtr.c([core_id],'vorticity_magnitude')


            avx = thtr.c([core_id],'velocity_x')
            avy = thtr.c([core_id],'velocity_y')
            avz = thtr.c([core_id],'velocity_z')
            abx = thtr.c([core_id],'magnetic_field_x')
            aby = thtr.c([core_id],'magnetic_field_y')
            abz = thtr.c([core_id],'magnetic_field_z')
            av2 = np.sqrt(avx**2 + avy**2 + avz**2)
            ab2 = np.sqrt(abx**2 + aby**2 + abz**2)
            vdotb = avx*abx + avy*aby + avz*abz
            costheta = vdotb/(av2*ab2)
            this_r = ms.r
            this_r[ this_r < min_r ] = min_r    
            for nf,frame in enumerate(thtr.frames):
                mask2 = ms.compute_unique_mask(core_id, dx=1./2048,frame=nf)
                density = all_density[mask2,nf]
                phi = all_phi[mask2,nf]
                cell_volume = all_cell_volume[mask2,nf]
                magnetic_field_strength = all_magnetic_field_strength[mask2,nf]
                kinetic_energy = all_kinetic_energy[mask2,nf]
                vorticity = all_vorticity[mask2,nf]

                r = this_r[mask2,nf]
                mass = (density*cell_volume).sum()
                volume = (cell_volume).sum()

                this_ctheta = costheta[mask2,nf]
                ctheta = (this_ctheta*cell_volume).sum()/volume


                radius_rms = np.sqrt((density*cell_volume*r**2).sum()/mass)
                #if Rhull is not None:
                #    radius_rms = Rhull
                G2 = -G*mass**2/radius_rms



                Eg = (density*phi*cell_volume).sum()
                Eb = (magnetic_field_strength*cell_volume).sum()
                Ek = (kinetic_energy*cell_volume).sum()
                W  = (vorticity*cell_volume).sum()
                self.Energy_G_core[core_id].append(Eg)
                self.Energy_B_core[core_id].append(Eb)
                self.Energy_M_core[core_id].append(G2)
                self.Energy_K_core[core_id].append(Ek)
                self.Energy_W_core[core_id].append(W)

                self.Energy_G_time[nf].append(Eg)
                self.Energy_B_time[nf].append(Eb)
                self.Energy_M_time[nf].append(G2)

                self.costheta[core_id].append(ctheta)

import three_loopers_energy as tl_energy
import convex_hull_tools as CHT
reload(CHT)
#if 'ht1' not in dir() or clobber: 
#    ht1 = CHT.hull_tool(tl.looper1)
if 'ht2' not in dir() or clobber:
    ht2 = CHT.hull_tool(tl_energy.looper2)
    ht2.make_hulls()
#if 'ht3' not in dir() or clobber:
#    ht3 = CHT.hull_tool(tl.looper3)

#import three_loopers as tl_full
if 'clobber' not in dir():
    clobber=True
#if 'mphi_tool1' not in dir() or clobber:
#    mphi_tool1=mphi_tool(tl_full.looper1)
#    mphi_tool1.run()
if 'mphi_tool2' not in dir() or clobber:
    mphi_tool2=mphi_tool(tl_energy.looper2)
    mphi_tool2.run(hull=ht2)
#if 'mphi_tool3' not in dir() or clobber:
#    mphi_tool3=mphi_tool(tl_full.looper3)
#    mphi_tool3.run()


if 0:
    #
    # rho Phi vs. gM^2/R, individual frames
    #
    fig,ax=plt.subplots(1,1, figsize=(12,4))
    axes=[ax]
    for nt,tool in enumerate([mphi_tool2]):# [mphi_tool1,mphi_tool2,mphi_tool3]):
        tff_global = np.sqrt(3*np.pi/(32*G*1))
        times = tool.this_looper.tr.times/tff_global
        for ntime, time in enumerate(times):
            ax.clear()
            the_x,the_y = np.log10(np.abs(tool.Energy_G_time[ntime])), np.log10(np.abs(tool.Energy_M_time[ntime]))
            ax.scatter( the_x, the_y)
            axbonk(ax, xlabel=r'$\int \rho \phi$', ylabel=r'$G M^2/R$', xlim=[-5,0], ylim=[-8,1])#, xscale='log',yscale='log')

            pfit = np.polyfit( the_x, the_y, 1)
            print(pfit)
            ax.plot( the_x, pfit[0]*the_x+pfit[1])
            #xline = nar([1e-3, 1e1])
            #yline = xline
            #ax.plot(xline,yline)
            outname = "plots_to_sort/Egrav_GMMr_n%04d.png"%ntime
            print(outname)
            fig.savefig(outname)

def derp(arr):
    return arr / np.abs(arr[:6].mean())


def load_quan(h5_name): 
    stuff={}
    if len(glob.glob(h5_name)): 
        fptr = h5py.File(h5_name,'r') 
        try: 
            for k in fptr: 
                stuff[k]=list(copy.copy(fptr[k].value)) 
        except: 
            raise 
        finally: 
            fptr.close() 
    else: 
        print("NO SUCH FILE", h5_name) 
    order = np.argsort( stuff['t'])
    for key in stuff:
        stuff[key] = nar(stuff[key])[order]

    return stuff

quan_2 = load_quan("/home/dccollins/ytscripts/quan_box_u202.h5")
#Mass to flux
#Field to Density
#Field to Vorticity
if 1:
    #
    # rho Phi gM^2/R, B^2 vs time
    #
    fig,ax=plt.subplots(2,4, figsize=(12,8))
    axes=ax[0]
    row1=ax[1]

    for nt,tool in enumerate([mphi_tool2]):# [mphi_tool1,mphi_tool2,mphi_tool3]):
        tff_global = np.sqrt(3*np.pi/(32*G*1))
        times = tool.this_looper.tr.times/tff_global
        name = tool.this_looper.out_prefix
        for ncore, core_id in enumerate(tool.cores_used):

            this_G = ((np.abs(tool.Energy_G_core[core_id])))
            this_B = ((np.abs(tool.Energy_B_core[core_id])))
            this_M = ((np.abs(tool.Energy_M_core[core_id])))
            this_K = ((np.abs(tool.Energy_K_core[core_id])))
            this_W = ((np.abs(tool.Energy_W_core[core_id])))


            axes[0].plot( times,derp(this_G))
            axes[1].plot( times,derp(this_B))
            axes[2].plot( times,derp(this_K))
            axes[3].plot( times,derp(this_W))
            axbonk(axes[0], xlabel=r'$t$', yscale='log', ylabel='Egrav')
            ylim = [2e-4,5]
            axbonk(axes[1], xlabel=r'$t$', yscale='log', ylabel='Emag', ylim=[2e-4,20])
            axbonk(axes[2], xlabel=r'$t$', yscale='log', ylabel='Ekin', ylim=ylim)
            axbonk(axes[3], xlabel=r'$t$', yscale='log', ylabel='Omega', ylim=ylim)



            row1[3].plot( quan_2['t']/tff_global, quan_2['Bfield_strength'],c='k')
            axes[1].plot( quan_2['t']/tff_global, quan_2['Bfield_strength']-quan_2['Bfield_strength'][-1],c='k')

            MassToFlux = np.sqrt( this_G/this_B)
            row1[0].plot(times, derp(MassToFlux))
            axbonk(row1[0], xlabel=r'$t$', yscale='log', ylabel=r'$M/\Phi$')

            FieldToVorticity = this_B/this_W
            row1[1].plot(times, FieldToVorticity)
            axbonk(row1[1], xlabel=r'$t$', yscale='log', ylabel=r'$B/\omega$')

            row1[2].plot(times, tool.costheta[core_id])
            #axbonk(row1[1], xlabel=r'$t$', yscale='log', ylabel=r'$B/\omega$')


        outname = "plots_to_sort/Egrav_GMMr_%s.png"%(name)
        print(outname)
        fig.savefig(outname)

if 0:
    fig,ax=plt.subplots(1,3, figsize=(12,4))
    axes=ax.flatten()
    for nt,tool in enumerate([mphi_tool2]):# [mphi_tool1,mphi_tool2,mphi_tool3]):
        G=1620./(4*np.pi)
        tff_global = np.sqrt(3*np.pi/(32*G*1))
        ntimes = tool.times.size
        ncores = len(tool.cores_used)
        masses = np.zeros([ntimes,ncores])
        these_times = tool.times/tff_global
        outname='plots_to_sort/mass_time_heatmap_unique_particles.pdf'
        for ncore,core_id in enumerate(tool.cores_used):
            this_mass = tool.unique_mass[core_id]
            if len(this_mass) == 0:
                pdb.set_trace()
            masses[:,ncore]=this_mass/nar(this_mass[:6]).mean()
            dof = nar(tool.dof[core_id])
            #masses[:,ncore] *= dof[0]/dof


        nc = len(tool.cores_used)
        take_a_few = ((nc-1)*np.random.random(10)).astype('int')
        print(take_a_few)
        for ncore,core_id in enumerate(nar(tool.cores_used)[take_a_few]):
            rel_mass=masses[:,ncore]

            axes[nt].plot(these_times, rel_mass,c=[0.5]*4)
            #axes[3].plot(these_times, tool.dof[core_id]/tool.dof[core_id][-1],c=[0.5]*4)

        mass_bins_edge = np.logspace(-3,4,101)
        mass_bins_edge = np.logspace(-2,1,101)
        xbins = these_times
        ybins = 0.5*(mass_bins_edge[1:]+mass_bins_edge[:-1])
        nx = len(xbins) ; ny=len(ybins)
        TheX = np.r_[(ny)*[xbins]].transpose()
        TheY = np.r_[(nx)*[ybins]]

        hist = np.zeros( [xbins.size,ybins.size])
        for ntime, time in enumerate(these_times):
            thishist,bins = np.histogram(masses[ntime,:],bins=mass_bins_edge)
            hist[ntime,:]=thishist

        cmap = copy.copy(mpl.cm.get_cmap("viridis"))
        cmap.set_under('w')
        minmin = hist[hist>0].min()
        norm = mpl.colors.LogNorm(vmin=1,vmax=33)
        ploot=axes[nt].pcolormesh(TheX, TheY, hist, cmap=cmap,norm=norm,shading='nearest')
        #axes[nt].plot(these_times, [2]*tool.times.size,c='r')#,lw=0.2)
        #axes[nt].plot(these_times, [1./2]*tool.times.size,c='r')#,lw=0.2)
        axbonk(axes[nt],ylabel=None,xlabel=r'$t/t_{\rm{ff}}$',yscale='log')
    axes[0].set_ylabel(r'$M/\bar{M_{\rm{early}}}$')
    fig.colorbar(ploot)
    fig.savefig(outname)
    print(outname)

