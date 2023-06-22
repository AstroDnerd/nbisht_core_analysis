from starter2 import *



this_simname = 't02'
frame = dl.target_frames[this_simname]
ds = yt.load("%s/DD%04d/data%04d"%(dl.sims[this_simname],frame,frame))

peak_fname = dl.peak_list[this_simname]
fptr = h5py.File(peak_fname,'r')
peaks = fptr['peaks'][()]
fptr.close()

axis = 0
proj = ds.proj('density',0)
pw = proj.to_pw()
for npeak,peak in enumerate(peaks):
    pw.annotate_text(peak,'.')
    pw.annotate_text(peak,r'$%d$'%npeak)
pw.save('plots_to_sort/%s_%04d_peak_loc'%(this_simname, frame))
