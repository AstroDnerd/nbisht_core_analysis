from starter2 import *
import sf2
reload(sf2)

if 'mysf' not in dir():
    mysf={}

sims = ['u201','u202','u203']
for sim in sims:
    if sim in mysf:
        continue

    frame = 0
    fname = dl.sims[sim] + "/DD%04d/data%04d"%(frame,frame)
    ds = yt.load(fname)

    mysf[sim]=sf2.make_sf(ds=ds)

fig,ax=plt.subplots(1,1)

color={'u05':'r','u10':'g','u11':'b'}
color={'u201':'r','u202':'g','u203':'b'}
for sim  in sims:
    rbins,SS = mysf[sim].bin_take3(); SS/=2*np.pi
    X = np.log10(rbins)
    Y = np.log10(SS)
    ok = Y>0
    ok = np.logical_and(ok, rbins<0.3)
    pfit = np.polyfit(X[ok],Y[ok],1)
    ax.plot(10**X[ok], 10**(X[ok]*pfit[0]+pfit[1]),c='k')

    Lsonic = 10**(-pfit[1]/pfit[0])

    ax.plot(rbins,SS,c=color[sim],label=r'$%s L_s = %0.4f$'%(sim,Lsonic))
#axbonk(ax,xlabel='r',ylabel='SF2',yscale='log',xscale='log')
ax.legend(loc=0)
fig.savefig('plots_to_sort/tmp.png')

