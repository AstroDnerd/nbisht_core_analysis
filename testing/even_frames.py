

plt.close('all')
fig,ax=plt.subplots(1,1)
times =  set_looper.tr.times
dtimes = times[1:]-times[:-1]
#ax.plot(dtimes,marker='*')
things = times % times[2]
print("N zeros", (things == 0).sum())
ax.hist( np.log10(things[things >0]), histtype='step')

fig.savefig('plots_to_sort/times.png')
