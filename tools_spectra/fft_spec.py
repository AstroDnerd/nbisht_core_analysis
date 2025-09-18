
frame =125
spec={}
directory="/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential"
fft = "%s/DD%04d.products/fft_density.float32"%(directory,frame)
thingy=dpy( fft, ['density'])[0]
density_cube=np.fft.ifft( thingy*np.conj(thingy))


