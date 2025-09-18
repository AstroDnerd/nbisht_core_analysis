from starter2 import *
import xtra_energy
plt.close('all')


frame=1
ddd='/data/cb1/Projects/P19_CoreSimulations/new_sims/u25b_little_big_2'
fff=ddd+'/DD%04d/data%04d'%(frame,frame)

ds = yt.load(fff)
xtra_energy.add_energies(ds)
cg = ds.covering_grid(0,[0.0,0.0,0.0],[128,128,128],fields = ['density', YT_grav_energy])

if 0:
    #projections
    fig,ax=plt.subplots(1,1)
    ax.imshow( cg[YT_density].sum(axis=0).v)
    fig.savefig('plots_to_sort/grav_20.png')

    proj=ds.proj(YT_grav_energy,0)
    proj.to_pw().save('plots_to_sort/grav_21.png')

if 0:
    cg.set_field_parameter('center',ds.arr([0.5]*3,'code_length'))
    fig,ax = plt.subplots(1,1)
    plt.scatter(cg['radius'],cg[YT_potential_field])
    fig.savefig('plots_to_sort/grav_22.png')

if 1:
    import gravity
    reload(gravity)
    G = ds.parameters['GravitationalConstant'] #wants the 4pi
    ggg = gravity.gravity(cg['density'].v, G)
    ggg.solve()

    fig,ax = plt.subplots(1,1)
    #ax.scatter(cg['radius'],ggg.phi)
    #ax.scatter(cg['radius'],cg[YT_potential_field])
    ax.scatter( ggg.phi, cg[YT_potential_field])

    fig.savefig('plots_to_sort/grav_23.png')

