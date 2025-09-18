

fig, ax = plt.subplot(2,2)
ax0=ax[0][0]; a1=ax[0][1]; a2=ax[1][0]; a3=ax[1][1]

#like np.array
x,y = np.mgrid[0:1:0.1, 0:1:0.2]
#like np.linspace.  Use a complex step size
x,y = np.mgrid[0:100:1j, 0:100:1j]


#For formatting latex tables:
a=np.arange(12)
"%s |"*a.size%tuple(a) 
