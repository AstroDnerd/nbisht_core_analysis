from starter2 import *

import gravity
reload(gravity)

import pcolormesh_helper as pch
reload(pch)

class cut_core():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.output = defaultdict(list)

    def run(self,core_list=None, do_plots=True,do_proj=True):
        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        this_looper=self.this_looper
        frame = this_looper.target_frame
        ds = this_looper.load(frame)
        G = ds['GravitationalConstant']/(4*np.pi)
        xtra_energy.add_energies(ds)
        for core_id in core_list:
            print('Potential %s %d'%(this_looper.sim_name,core_id))

            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            c = nar([ms.mean_x[-1], ms.mean_y[-1],ms.mean_z[-1]])
            
            R_SPHERE = 32/128
            left = c - R_SPHERE/2
            right = c+R_SPHERE/2
            dx =  1/2048
            size = ((right-left)/dx).astype('int')
            cg = ds.covering_grid(4, left, [512]*3, fields=[YT_density,YT_potential_field])
            
            print('do the solve')
            
            ggg = gravity.gravity( cg['density'], ds.parameters['GravitationalConstant'])
            ggg.solve()
            print('done')
            if 1:
                #plot density
                print('proj')
                fig,ax=plt.subplots(1,1)
                norm=mpl.colors.LogNorm(vmin=cg[YT_density].min(),vmax=cg[YT_density].max())
                ax.imshow( cg[YT_density].sum(axis=0),norm=norm)
                fig.savefig('plots_to_sort/grav_rho_%s_c%04d.png'%(self.this_looper.sim_name,core_id))

            if 1:
                #plot phi
                print('proj2')
                fig,ax=plt.subplots(1,2,figsize=(12,8))
                phi=cg[YT_potential_field].sum(axis=0)
                norm = mpl.colors.Normalize( vmin = phi.min(), vmax=phi.max())
                ploot=ax[0].imshow( phi,norm=norm)
                ax[0].set_title('phi disk')
                fig.colorbar(ploot,ax=ax[0])


                #my phi
                print('proj3')
                phi=ggg.phi.sum(axis=0)
                norm = mpl.colors.Normalize( vmin = phi.min(), vmax=phi.max())
                ploot=ax[1].imshow( phi,norm=norm)
                ax[1].set_title('phi me')
                fig.colorbar(ploot,ax=ax[1])

                fig.savefig('plots_to_sort/grav_phi_%s_c%04d.png'%(self.this_looper.sim_name,core_id))

            plt.close('all')
            fig,ax=plt.subplots(1,2)
            cg.set_field_parameter('center',ds.arr(c,'code_length'))
            rbins = np.geomspace( 1/2048, cg['radius'].max(),128)
            phibins=np.linspace( cg[YT_potential_field].min(), cg[YT_potential_field].max(),125)
            #Bbins = np.geomspace( Beta[Beta>0].min(), d.max(),64)
            #Lbins = np.geomspace( Lmag[Lmag>0].min(), b.max(),64)
            rrr=cg['radius'].v.flatten()
            ppp=cg[YT_potential_field].v.flatten()
            ddv=cg[YT_cell_volume].v.flatten()
            hist, xb, yb = np.histogram2d( rrr, ppp, bins=[rbins,phibins],weights=ddv)
            pch.helper(hist,xb,yb,ax=ax[0])

            phibins=np.linspace( ggg.phi.min(), ggg.phi.max(),128)
            ppp=ggg.phi.flatten()
            hist, xb, yb = np.histogram2d( rrr, ppp, bins=[rbins,phibins],weights=ddv)
            pch.helper(hist,xb,yb,ax=ax[1])
            ax[0].set_xscale('log')
            ax[1].set_xscale('log')

            #ax.scatter( cg['radius'],np.abs(cg[YT_potential_field]))
            #ax.scatter( cg['radius'],np.abs(ggg.phi))
            #ax.set_yscale('log')
            fig.savefig('plots_to_sort/phi_scatter_%s_c%04d.png'%(self.this_looper.sim_name,core_id))



            print(size)
            if 0:
                rsph = ds.arr(R_SPHERE,'code_length')
                sp = ds.sphere(c,rsph)
                GE = np.abs(sp[YT_grav_energy])
                dv = np.abs(sp[YT_cell_volume])
                RR = sp[YT_radius]
                DD = sp[YT_density]
            

import three_loopers_six as TL
simlist=['u601']
for sim in simlist:
    all_cores=np.unique( TL.loops[sim].tr.core_ids)
    #core_list=list(all_cores)
    #core_list=all_cores[:1]
    core_list=[323]
    obj = cut_core(TL.loops[sim])
    obj.run(core_list=core_list)


