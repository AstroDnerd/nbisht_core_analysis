
this_tr =TL.looper1.tr 
dx = 1./2048
ix =np.floor(this_tr.c([core_id],'x')/dx)
iy =np.floor(this_tr.c([core_id],'y')/dx)
iz =np.floor(this_tr.c([core_id],'z')/dx)
index = ix + nx*(iy * nx*iz)
ar = np.argsort(index,axis=0)
rs = np.argsort(ar,axis=0)
isorted=np.take_along_axis(index,ar,axis=0)
mask = np.ones_like(ix,dtype='bool')
mask[1:,:] = isorted[1:,:]-isorted[:-1,:] != 0
mask2 = np.take_along_axis(mask, rs, axis=0)
