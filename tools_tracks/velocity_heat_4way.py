from starter2 import *

import hair_dryer
reload(hair_dryer)


import three_loopers_u500 as tl5
simlist = ['u501','u502','u503']
simlist = ['u503']
loops = tl5.loops


for this_simname in simlist:
    this_looper=loops[this_simname]
    core_list=np.unique(this_looper.tr.core_ids)


    hair = hair_dryer.hair_tool(loops[this_simname])

    for core_id in core_list:
        fig,ax=plt.subplots(2,2,figsize=(12,10))
        axlist = ax.flatten()

        hair.run(core_list=[core_id],external_ax=ax[0][0])
        tff=r'$t/t_{\rm{ff}}$'
        axbonk(ax[0][0],xlabel=r'y [code_units]', ylabel='z [code_units]')

        for ns, field in enumerate(['velocity_x','velocity_y','velocity_z']):
            this_ax = axlist[ns+1]
            bins = np.linspace(-16,16,65)
            heat_map.heat_for_quantity( this_looper, field=field,core_list=[core_id],external_ax=this_ax, bins=bins)
            this_ax.plot([0.8,0.8],[12,13],c='r',linewidth=5)
            axbonk(this_ax,xlabel=tff,ylabel=field)

        outname='plots_to_sort/hair_velocity_4way_%s_c%04d.png'%(this_simname,core_id)
        fig.savefig(outname)
        print(outname)


        plt.close(fig)


