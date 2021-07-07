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

subset = 2
for sim in sims:
    if sim in msf:
        continue
    directory = dl.sims[sim]
    ds = yt.load("%s/DD%04d/data%04d"%(directory,frame,frame))
    cg = ds.covering_grid(0,ds.domain_left_edge,ds.domain_dimensions)
    msf[sim] = sf2l.get_sf2l(cg, subset=subset, name=sim)


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
    binnedplot(fig,ax[ns],Rcube, SF2cube)
    save(Rcube,SF2cube,'SF2_L_subset_%d_%s.h5'%(subset,sim))

fig.savefig('plots_to_sort/all_sf2L.png')


