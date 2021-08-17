


fig,ax=plt.subplots(1,1)


frame_ind = np.where(this_looper.tr.frames == frame)[0][0]
density3 = this_looper.tr.c(core_list,FIELD)[:,frame_ind]
cell_volume3 = this_looper.tr.c(core_list,'cell_volume')[:,frame_ind]
#density3 = this_looper.tr.track_dict[FIELD][:,frame_ind]
#cell_volume3 = this_looper.tr.track_dict['cell_volume'][:,frame_ind]
vals3, bbb3 = np.histogram(density3, weights=cell_volume3, bins=bbb1)
bcen3=0.5*(bbb3[1:]+bbb3[:-1])
db3 = bbb3[1:]-bbb3[:-1]

mask = (ad[deposit_tuple] > 0 )
density4 =  ad['velocity_magnitude'][mask] 
cell_volume4 = ad['cell_volume'][mask]

vals4, bbb4 = np.histogram(density4, weights=cell_volume4, bins=bbb1, density=True)
bcen4=0.5*(bbb4[1:]+bbb4[:-1])
db4 = bbb4[1:]-bbb4[:-1]


ax.plot( bcen1, vals1,'k')  #the full hist
ax.plot( bcen2, vals2,'r')  #what yt givs
#ax.plot( bcen3, vals3,'g')  #from trackage
ax.plot( bcen3, vals3/db3,'b')  #from trackage
#ax.plot( bcen4, eta1*vals4,c='k')
axbonk(ax,xscale='log',yscale='log')
fig.savefig('plots_to_sort/vel_test.png')

#mask = (ad[deposit_tuple] > 0 )
#v_m = np.sort( ad['velocity_magnitude'][mask] )
#v_t_2 = np.sort( this_looper.tr.track_dict['velocity_magnitude'][:,0])
##v_t = np.sort( this_looper.tr.c(core_list,'velocity_magnitude')[:,0] )
