from starter2 import *
#import three_loopers_1tff as tl
import sf2l 
reload(sf2l)
import radial_binner as rb

fig,ax=plt.subplots(1,3, figsize=(12,4))
frame=0
sims=['u05','u10','u11']
if 'msf' not in dir():
    msf = {}

if 'vrms' not in dir():
    vrms = {}


for sim in sims:
    if sim in vrms:
        continue
    directory = dl.sims[sim]
    ds = yt.load("%s/DD%04d/data%04d"%(directory,frame,frame))
    cg = ds.covering_grid(0,ds.domain_left_edge,ds.domain_dimensions)
    print('vrms')
    vx = cg['velocity_x']
    vy = cg['velocity_y']
    vz = cg['velocity_z']
    vx0 = vx.mean()
    vy0 = vy.mean()
    vz0 = vz.mean()

    vrms[sim] = np.sqrt( (vx-vx0)**2+(vy-vy0)**2+(vz-vz0)**2)

for sim in sims:
    if sim in msf:
        continue
    directory = dl.sims[sim]
    ds = yt.load("%s/DD%04d/data%04d"%(directory,frame,frame))
    cg = ds.covering_grid(0,ds.domain_left_edge,ds.domain_dimensions)
    msf[sim] = sf2l.get_sf2l(cg, subset=0)

def image(fig,ax,Rmag,SF2):

    #x0 = -0.5
    #dx=1./128
    #Rmag2=np.mgrid[0.5*dx+x0:1+x0:dx,0.5*dx+x0:1+x0:dx,0.5*dx+x0:1+x0:dx]
    #RmagN=  np.sqrt(Rmag2[0,...]**2+Rmag2[1,...]**2+Rmag2[2,...]**2)
    ok = SF2>0
    Xfield = Rmag[ok].flatten()/128
    Yfield = SF2[ok].flatten()

    rbins = np.linspace(0,1,64)
    minmin = Yfield[ Yfield>0].min()
    Sbins = np.logspace( np.log10(minmin), np.log10(Yfield.max()))
    bins=[rbins,Sbins]

    this_hist, xedge, yedge= np.histogram2d(Xfield,Yfield)#,bins=bins)

    minmin_hist = this_hist[ this_hist>0].min()
    norm = mpl.colors.LogNorm(vmin=minmin_hist,vmax=this_hist.max())
    nx = xedge.size-1
    ny = yedge.size-1
    def cen(arr):
        return 0.5*(arr[1:]+arr[:-1])
    TheX = np.r_[(ny)*[cen(xedge)]].transpose()
    TheY = np.r_[(nx)*[cen(yedge)]]
    cmap = copy.copy(mpl.cm.get_cmap("viridis"))
    cmap.set_under('w')
    ax.pcolormesh(TheX, TheY, this_hist, cmap=cmap, norm=norm,shading='nearest')

def binnedplot(fig,ax,Rmag,SF2):
    dx = 1./128
    bins = np.arange(0.5*dx,2,dx)
#reload(rb)
    binner = rb.rb2(Rmag/128, SF2, bins=bins)
    rbins = binner[1]
    SF2bins = binner[2]
    ax.plot( rbins, SF2bins)

    ok = np.logical_and(rbins>0, SF2bins>0)
    ok = np.logical_and(ok, np.isnan(SF2bins)==False)
    ok = np.logical_and(ok, rbins < 0.5)
    
    if 0:
        X, Y = binner[1][ok], binner[2][ok]
        pfit = np.polyfit(np.log10(X),np.log10(Y),1)
        logx = np.log10(X)
        ax.plot( X, 10**(logx*pfit[0]+pfit[1]),c='k')
    if 1:
        X, Y = binner[1][ok], binner[2][ok]
        pfit = np.polyfit(X,Y,1)
        logx = np.log10(X)
        ax.plot( X, (X*pfit[0]+pfit[1]),c='k')



def dump_SF2(Rmag, SF2, outname):

    dx = 1./128
    bins = np.arange(0.5*dx,2,dx)
    binner = rb.rb2(Rmag/128, SF2, bins=bins)
    ax.plot( binner[1], binner[2])
    dpy_save(outname,{'R':binner[1],'SF2L':binner[2]},fields=['R','SF2L'])
    fig.savefig('plots_to_sort/thing2.png')
    plt.close(fig)

def save(Rmag, SF2, outname):
    fptr=h5py.File(outname,'w')
    fptr.create_dataset('Rmag',data=Rmag)
    fptr.create_dataset('SF2_L',data=SF2)
    fptr.close()
fig,ax=plt.subplots(1,3,figsize=(12,4))

for ns,sim in enumerate(sims):
    Rcube, SF2cube = msf[sim][0],msf[sim][1]
    SF2cubenorm = SF2cube/SF2cube.size
    plot(fig,ax[ns],Rcube,SF2cube)
    binnedplot(fig,ax[ns],Rcube, SF2cube)
    save(Rcube,SF2cube,'SF2_L_%s.h5'%sim)

fig.savefig('plots_to_sort/all_sf2L.png')


