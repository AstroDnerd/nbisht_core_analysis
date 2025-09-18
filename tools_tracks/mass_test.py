

rng = 10
Np = 100
index = np.floor(rng*np.random.random(Np) )
density = 10*index

ar = np.argsort(index)
rs = np.argsort(ar)
isorted=index[ar]
mask = np.ones_like(density,dtype='bool')
mask[1:] = isorted[1:]-isorted[:-1] != 0
mask2 = mask[ rs]
mass = (density[mask2])

