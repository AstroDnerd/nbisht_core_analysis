
from starter2 import *

plot_dir = './plots_to_sort'
for sim_name in ['u201','u202','u203']:
    directory = dl.sims[sim_name]
    frame = dl.target_frames[sim_name]
    set_name = "%s/DD%04d/data%04d"%(directory,frame,frame)
    ds = yt.load(set_name)
    proj = ds.proj('density',0)
    pw = proj.to_pw(width = (1.0,'code_length'), origin='domain')
    pw.set_cmap('density','Greys')
    pw.set_axes_unit('code_length')
    pw.save('%s/%s_proj'%(plot_dir,sim_name))

