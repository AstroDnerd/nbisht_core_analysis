
import os
import numpy as na
def read_alexei(in_file, Isothermal=False,dtype='<f4',dtype_out='float32'):
    """Reads a set of files as given by Alexei.
    Converts pressure to energy"""
    mv = 8
    if Isothermal:
        mv = 7
    gamma = 5./3    
    d_file  = open(in_file+'.dat1','rb')
    vx_file = open(in_file+'.dat2','rb')
    vy_file = open(in_file+'.dat3','rb')
    vz_file = open(in_file+'.dat4','rb')
    bx_file = open(in_file+'.dat5','rb')
    by_file = open(in_file+'.dat6','rb')
    bz_file = open(in_file+'.dat7','rb')
    if Isothermal == False:
        p_file  = open(in_file+'.dat8','rb')
    else:
        p_file = None

    count_debug = -1
    try:
       d  = na.fromfile(d_file ,dtype=dtype,count = count_debug) 
       print "read density"
       vx = na.fromfile(vx_file,dtype=dtype,count = count_debug) 
       print "read vx"
       vy = na.fromfile(vy_file,dtype=dtype,count = count_debug) 
       print "read vy"
       vz = na.fromfile(vz_file,dtype=dtype,count = count_debug) 
       print "read vz"
       bx = na.fromfile(bx_file,dtype=dtype,count = count_debug) 
       print "read bx"
       by = na.fromfile(by_file,dtype=dtype,count = count_debug) 
       print "read by"
       bz = na.fromfile(bz_file,dtype=dtype,count = count_debug) 
       print "read bz"
       if Isothermal == False:
           p  = na.fromfile(p_file ,dtype=dtype,count = count_debug) 
           e =  p/(gamma-1) + 0.5*d*(vx*vx+vy*vy+vz*vz) + 0.5*(bx*bx+by*by+bz*bz)
           print "read p"
       else:
           p = 0
           e = 0

    except:
        raise
    finally:
        d_file.close()
        vx_file.close()
        vy_file.close()
        vz_file.close()
        bx_file.close()
        by_file.close()
        bz_file.close()
        if Isothermal == False:
            p_file.close()

    #pdb.set_trace()
    if Isothermal == False:
        field_list = [d,vx,vy,vz,bx,by,bz,e,p]
    else:
        field_list = [d,vx,vy,vz,bx,by,bz]
    return field_list
