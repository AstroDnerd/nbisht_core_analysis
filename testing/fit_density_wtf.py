from starter2 import *
frame=0
plt.close('all')
#xtra_energy.add_energies(ds)
#xtra_energy.add_gdotgradrho(ds)
#nf = np.where( this_looper.tr.frames == frame)[0][0]
#time = times[nframe,0]

#c = nar([ms.mean_xc[nf], ms.mean_yc[nf],ms.mean_zc[nf]])
#msR = ms.rc
#msR[ msR<1/2048]=1/2048
#
#MaxRadius=msR[:,nf].max()
#Radius = max([8.0/128, MaxRadius])

fig6,ax6=plt.subplots(3,2)
for NROW in [0]:
    c1=nar([0.40912891, 0.15447959, 0.0391341 ])
    rsph=0.0625
    c = np.round(c1*4096)/4096
    print((c-c1)*2048)

    if NROW==-1:
        fname='/data/cb1/Projects/P19_CoreSimulations/u202-Beta2/GravPotential/DD0090/data0090'
        #c=[0.40912891, 0.15447959, 0.0391341 ]
        #rsph=0.0625/4
        ds=yt.load(fname)
        sp = ds.sphere(c,rsph)
        dv = sp[YT_cell_volume]
        RR = sp['radius']
        DD = sp[YT_density]
        if 1:
            proj=ds.proj(YT_density,0,data_source=sp, center=c)
            pw=proj.to_pw(center=c)
            pw.zoom(4)
            pw.save('plots_to_sort/cg1')
    if NROW==0:
        fname='/data/cb1/Projects/P19_CoreSimulations/u202-Beta2/GravPotential/DD0090/data0090'
        L = c-rsph
        R = c+rsph
        dsreg=yt.load(fname)
        reg = dsreg.region(c,L,R)
        dv = reg[YT_cell_volume]
        RR = reg['radius']
        DD = reg[YT_density]
        if 1:
            proj=ds.proj(YT_density,0,data_source=reg, center=c)
            pw=proj.to_pw(center=c)
            pw.zoom(1/(R[0]-L[0]))
            pw.save('plots_to_sort/cg2')

    if NROW==1:
        fname='/data/cb1/Projects/P19_CoreSimulations/u202-Beta2/GravPotential/DD0090/data0090'
        if 0:
            rsph=0.0625/16
            L = c-rsph
            R = c+rsph
            #L =nar([0.34662891, 0.09197959, 0.9766341 ])
            #R =nar([0.47162891, 0.21697959, 1.1016341 ])
        L =reg.left_edge.v-0.5/2048
        R =reg.right_edge.v-0.5/2048
        ds=yt.load(fname)
        #sp = ds.sphere(c,rsph)
        level=4
        dx = ds.index.get_smallest_dx()
        Nzones = ((R-L)/dx).astype('int')
        cg = ds.covering_grid(level,L,Nzones)
        fig0,ax0=plt.subplots(1,1)
        ax0.imshow( np.log10( cg['density'].sum(axis=0)).transpose(), origin='lower')
        fig0.savefig('plots_to_sort/cg0')

        


        dv = cg[YT_cell_volume].flatten()
        RR = np.sqrt( (cg[YT_x].v-c[0])**2+(cg[YT_y].v-c[1])**2+(cg[YT_z].v-c[2])**2)   
        ok = RR>0
        DD = cg[YT_density][ok].flatten()
        fig1,ax1=plt.subplots(1,1)
        ax1.imshow( RR.sum(axis=0))
        fig1.savefig('plots_to_sort/fart')
        RR = RR[ok].flatten()


    if NROW==2:
        x,y,z=np.mgrid[-0.5:0.5:1/16, -0.5:0.5:1/16, -0.5:0.5:1/16]
        RR = np.sqrt(x**2+y**2+z**2).flatten()
        DD = np.zeros_like(RR)
        DD[RR>0] = RR[RR>0]**-1.39
        dv = np.zeros_like(RR)+1/16**3
        ok = RR>0
        RR = RR[ok]
        DD = DD[ok]


    figa,axa=plt.subplots(1,1)
    thing1=nar(sorted( cg[YT_x].flatten().v))
    thing2=nar(sorted( reg[YT_x].v))
    axa.plot( thing1[:10],c='g')
    axa.plot( thing2[:10],c='r')
    figa.savefig('plots_to_sort/dump')



    ax1=ax6[NROW][0]
    ax2=ax6[NROW][1]
    ax1.set(title='rho')
    ax2.set(title='mass')


    ORDER = np.argsort( RR)
    RR_cuml = RR[ORDER]
    rho_srt = DD[ORDER]
    M_cuml = np.cumsum( DD[ORDER]*dv[ORDER])
    dv_sort = dv[ORDER]

    r_bins=np.geomspace(RR.min(),RR.max(),64)
    r_cen=0.5*(r_bins[1:]-r_bins[:-1])
    digitized = np.digitize( RR_cuml, r_bins)
    def digit(arr):
        output  =nar([ arr[ digitized == i].mean() if (digitized==i).any() else np.nan for i in range(1,len(r_bins))])
        return output


    pfit_rho = np.polyfit( np.log10(RR_cuml), np.log10(rho_srt),1, w=dv_sort)
    ax1.plot( RR_cuml, rho_srt)
    ax1.plot( RR_cuml, 10**(pfit_rho[0]*np.log10(RR_cuml)+pfit_rho[1]),label='a %0.2f'%(pfit_rho[0]))
    ax1.legend(loc=0)
    ax1.set(title='rho %0.2f'%(pfit_rho[0]))

    rho_bin = digit(rho_srt)
    ax1.plot( r_cen, rho_bin,c='k')
    #ax1.set(xscale='log',yscale='log')

    pfit_mass = np.polyfit( np.log10(RR_cuml[1:]), np.log10(M_cuml[1:]),1)

    MMM = pfit_rho[0]+3
    BBB = pfit_mass[1]
    ax2.plot( RR_cuml, M_cuml)
    ax2.plot( RR_cuml, 10**(pfit_mass[0]*np.log10(RR_cuml)+pfit_mass[1]),label='fit %0.2f'%( pfit_mass[0]))
    ax2.plot( RR_cuml, 10**(MMM*np.log10(RR_cuml)+BBB),label='a+3: %0.2f'%(MMM), c='r')
    ax2.set(title='M')
    ax2.legend(loc=0)

#ax6[1][0].plot(RR_cuml, rho_toy*4*np.pi*RR_cuml**3)
    for aaa in ax6.flatten():
        aaa.set(yscale='log',xscale='log')

fig6.savefig('plots_to_sort/derp3.png')
