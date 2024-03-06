"""
code for historical purposes only.
Don't run or import.
"""


FAIL NOW.
import os
import numpy as na
def read_alexei(in_file, out_dir,out_prefix,read_only=False,refine_by=0,
               Isothermal=False,dtype='<f4',dtype_out='float32'):
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
    if read_only:
        return field_list
        
    if 1:
        filename_dat = "%s/%s.dat"%(out_dir,out_prefix)
        file=open(filename_dat,'wb')
        print "write file",filename_dat
        field_list_2=[]
        for i,field in enumerate(field_list[0:8]):
            if refine_by:
                N = int( field.size**(1./3) +0.5 )
                field.shape = (N,N,N)
                print "Sub Sample: ", field.shape, "by",refine_by
                if i == 0:
                    SmallDensity = volavg(field,3,refine_by)
                    out_field = SmallDensity
                elif i in [1,2,3]:
                    out_field = volavg(d*field,3,refine_by)/SmallDensity
                else:
                    out_field  = volavg(field,3,refine_by)
                stat(out_field,'Writing %d'%(i))
            else:
                out_field = field

            field_list_2.append(out_field)
            if i < 7:
                """we write pressure."""
                file.write( out_field.astype(dtype_out))
            print 'write field',i,'min =',field.min()
        if Isothermal == False:
            d,u,v,w,a,b,c,t=field_list_2
            p=(gamma-1)*(t - 0.5*d*(u*u+v*v+w*w) - 0.5*(a*a+b*b+c*c) )
            file.write( p.astype(dtype_out) )
            
        file.close()
        #pdb.set_trace()
        filename_g2g = "%s/%s.g2g"%(out_dir,out_prefix)
        file = open(filename_g2g,'w')
        print filename_g2g
        file.write(" &in  file='%s.dat' offset=0.0 do_exp=f /\n"%(out_prefix))
        N = int( field.size**(1./3) +0.5 )
        if refine_by:
            size = nint(1.*N/refine_by)
        else:
            size = N
        file.write(" &out file='%s_regrid.dat' m=%d mv=%d offset=-.5 do_smooth=f /\n"%(out_prefix,size,mv))
        file.close()
          
    if 0:
        field_names = ['Density','Vx','Vy','Vz','Bx','By','Bz','E']
        out_filenames = ['%s/%s_%s'%(out_dir,out_prefix,name) for name in field_names] 
        for i in range(len(field_list)):
            print "Write field",field_names[i]
            out_grid = h5py.File(out_filenames[i],'w')
            try:
                F = field[i].swapaxes(0,2)
                out_grid.create_dataset(out_prefix+"_"+field_names[i],F.shape,data=F)
            except:
                raise
            finally:
                out_grid.close()

def aake_binary(oober,LiCode=False):
    if LiCode:
        fields=['Density','x-velocity','y-velocity','z-velocity',
                'MagneticField_C_1','MagneticField_C_2','MagneticField_C_3']
    else:
        fields=['Density','x-velocity','y-velocity','z-velocity','Bx','By','Bz']
    max_level=0
    for n in oober.frames:
        oober.fill(n)
        h=oober.h
        right=na.amax(h.gridRightEdge,axis=0)
        left=na.amin(h.gridLeftEdge,axis=0)
        print left, right
        resolution = (right-left)/h.grids[0]['dx']
        #cg=h.covering_grid(0,left,right,resolution,fields=fields)
        cg=h.covering_grid(0,left,resolution,fields=fields)

        return cg
def write_binary(oober,fields=None,LiCode=False,Alexei=False,Aake=True,AakeMeta=True,AakePower=True,cg_input=None,
                num_ghost_zones=0,format='float64'):
    """Write binary files from enzo data.
    Writes one for each *n* in *oober*.frames
    *LiCode* uses naming convention from the Enzo MHDCT version
    *Alexei* makes 7 data cubes, one for each field.
    *Aake* make one data cube with all of them.
    """
    if fields != None:
        fields = ensure_list(fields)
    if fields != None and len(fields) == 1:
        OneField = True
    else:
        OneField = False
    if fields == None:
        if LiCode:
            fields=['Density','x-velocity','y-velocity','z-velocity',
                    'MagneticField_C_1','MagneticField_C_2','MagneticField_C_3']
        else:
            fields=['Density','x-velocity','y-velocity','z-velocity','Bx','By','Bz']
    max_level=0
    for n in oober.frames:
        oober.fill(n)
        directory=oober.directory+"/DD%04d.products"%(n)
        if glob.glob(directory) == []:
            print "making directory",directory
            os.mkdir(directory)
        h=oober.h
        right=na.amax(h.grid_right_edge,axis=0)
        left=na.amin(h.grid_left_edge,axis=0)
        print left, right
        resolution = (right-left)/h.grids[0]['dx']+2*num_ghost_zones
        #cg=h.covering_grid(0,left,right,resolution,fields=fields)
        if (Aake or Alexei) and cg_input == None:
            cg=h.covering_grid(0,left,resolution,num_ghost_zones=num_ghost_zones)#,fields=fields)
        if cg_input != None:
            cg = cg_input
        if Aake:
            if OneField:
                filename = oober.directory+"/DD%04d.products/data%04d_%s.dat"%(n,n,fields[0])
            else:
                filename = oober.directory+"/DD%04d.products/data%04d.dat"%(n,n)
            file=open(filename,'wb')
            print "write file",filename
            for field in fields:
                file.write( cg[field].astype(format) )
            file.close()
        if AakeMeta:
            if OneField:
                filename = oober.directory+"/DD%04d.products/data%04d_%s.dim"%(n,n,fields[0])
            else:
                filename = oober.directory+"/DD%04d.products/data%04d.dim"%(n,n)
            file = open(filename,'w')
            file.write('&DIM\n')
            file.write('MX=%d, MY=%d, MZ=%d,'%tuple(resolution)) 
            file.write('MV=7, SX=1.00, SY=1.00, SZ=1.0, OFFSET=-0.50, X_INDEX=2, 4, DO_EXP=F, DO_MASSFLUX=F/')
            file.close()
            print filename
        if AakePower:
            if OneField:
                filename =  oober.directory+"/DD%04d.products/g2p_%s.in"%(n,fields[0])
                datafile = "data%04d_%s.dat"%(n,fields[0])
                outfile = 'power_%s.out'%fields[0]
            else:
                filename =  oober.directory+"/DD%04d.products/g2p.in"%n
                datafile = "data%04d.dat"%(n)
                outfile = 'power.out'
            file=open(filename,'w')
            file.write("&in\n infile='%s'\n outfile='%s'\n compensate=0.\n  do_average=t\n /"%(datafile,outfile))
            file.close()
                    
        if Alexei:
            for m,field in enumerate(fields):
                filename = directory + "/cube%04d.dat%d"%(n,m+1)
                print filename, "field",field
                file=open(filename,'wb')
                file.write( cg[field].astype('float32') )
                file.close()
    
    return cg
if 0:
    kitp_256_li = uber(dir='/scratch/david_c/Paper10/j14_256_rage/ic_take2',frames=[0],name='kitp_256_li',OutputName='MyOutputs.FileStaticOutputRaw')
    write_binary(kitp_256_li,LiCode=True,format='float32')
if 1:
    kitp_256_li = uber(dir='/scratch/david_c/Paper10/Enzo256_Hancock',
                       frames=[0,1],name='kitp_256_hancock',OutputName='load')
    write_binary(kitp_256_li,LiCode=True,format='float32')

if 0:
    kitp_512_li = uber(dir='/scratch/david_c/Paper10/j16_512',frames=[0],name='kitp_512_li',OutputName='MyOutputs.FileStaticOutputRaw')
    write_binary(kitp_512_li,LiCode=True,format='float32')
