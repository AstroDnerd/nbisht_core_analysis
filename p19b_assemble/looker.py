from starter1 import *


f1 = '/anvil/scratch/x-ux454321/p78c_high_res/IC_assembler/Turb_64_ic/kitp_ath_64.Bx'
f2 = '/anvil/scratch/x-ux454321/p78c_high_res/IC_assembler/RUN32/bx'
fnames = [f1,f2]

fptrs=[h5py.File(fn) for fn in fnames]
