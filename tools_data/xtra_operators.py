from starter2 import *
def gradf(field,direction,dds):
    iM1 = slice(None,-2)
    iP1 = slice(2,None)
    all = slice(1,-1)
    all_all=[all]*3
    dxi=1./(2*dds )
    out = np.zeros_like(field*dxi[0])
    Left = [all]*3
    Right = [all]*3
    Right[direction] = iP1
    Left[direction] = iM1
    Left = tuple(Left); Right=tuple(Right);all_all=tuple(all_all)
    out[all_all] = (field[ Right ]- field[ Left]) *dxi[direction]
    return out
def gradf2d(field,direction,dds):
    iM1 = slice(None,-2)
    iP1 = slice(2,None)
    all = slice(1,-1)
    all_all=[all]*2
    Left = [all]*2
    Right = [all]*2
    dxi=1./(2*dds )
    out = np.zeros_like(field*dxi[0])
    Right[direction] = iP1
    Left[direction] = iM1
    Left = tuple(Left); Right=tuple(Right);all_all=tuple(all_all)
    out[all_all] = (field[ Right ]- field[ Left]) *dxi[direction]
    return out
def grad(data,fieldname,direction):
    iM1 = slice(None,-2)
    iP1 = slice(2,None)
    all = slice(1,-1)
    all_all=tuple([all]*3)
    dxi=1./(2*data.dds )
    out = np.zeros_like(data[fieldname]*dxi[0])
    Left = [all]*3
    Right = [all]*3
    Right[direction] = iP1
    Left[direction] = iM1
    Left = tuple(Left); Right=tuple(Right)
    out[all_all] = (data[fieldname][ Right ]- data[fieldname][ Left]) *dxi[direction]

    return out

def AdotDel(data,fields1,field2):
    """returns A \cdot \nabla B
    where A (fields1) is a vector field and B (fields2) is a scalar field.
    """
    iM1 = slice(None,-2)
    iP1 = slice(2,None)
    all = slice(1,-1)
    all_all=[all]*3
    dxi=1./(2*data.dds )
    out = np.zeros_like(data[field2]*data[fields1[0]]*dxi[0])
    temp = np.zeros_like(data[field2])
    for i, fi in enumerate(fields1):
        Left = [all]*3
        Right = [all]*3
        Right[i] = iP1
        Left[i] = iM1
        Left = tuple(Left); Right=tuple(Right);all_all=tuple(all_all)
        temp[all_all] = data[field2][ Right ]- data[field2][ Left] 
        out[all_all] += data[fi][all_all]*(temp[all_all])*dxi[i]
    return out
    
def tensor_derivative(data,fields1,fields2):
    """performs (A \cdot \nabla) B,
    specifically, 
    Ai di Bj
    where di = d/d[x,y,z]
    Ai in fields1, Bj in fields2
    Assumes centered difference.
    !!!! I think there was an error, need to included - in out[j][all_all] line.  has been added, please check."""
    iM1 = slice(None,-2)
    iP1 = slice(2,None)
    all = slice(1,-1)
    all_all=[all]*3
    out = [np.zeros(data[fields1[0]].shape)]*3
    dxi=1./(2*data.dds )
    for j,fj in enumerate(fields2):
        for i, fi in enumerate(fields1):
            Left,Right = [all]*3,[all]*3
            Right[i] = iP1
            Left[i] = iM1
            out[j][all_all] += data[fi][all_all]*(data[fj][ Right ]- data[fj][ Left] )*dxi[i]
    return out

def GradScalar(V, dx_in):
    """takes the gradiant of a scalar"""
    #slice(start,stop,step), None defaults to everything
    if True:
        #Centered difference: dv/dx_{i} = 1/(2 dx)(v_{i+1} - v_{i-1})
        SL = slice(None,-2,None)  #V_{i-1}
        SR = slice(2,None,None)   #V_{i+1}
        act=slice(1,-1,None)      #active region
        width = 2
    dx = nar(dx_in)*width
    output = [np.zeros( V.shape ),np.zeros( V.shape ),np.zeros( V.shape ) ]
 
    #dV/dx \hat{x}
    dVdX = 1./dx[0]*(V[SR,act,act] - V[SL,act,act])
    output[0][act,act,act] = dVdX

    #dV/dy \hat{y}
    dVdX = 1./dx[1]*(V[act,SR,act] - V[act,SL,act])
    output[1][act,act,act] = dVdX

    #dV/dz \hat{z}
    dVdX = 1./dx[2]*(V[act,act,SR] - V[act,act,SL])
    output[2][act,act,act] = dVdX
    return output


def Curl(V,dx_in,component=None):
    """takes the curl of vector field V.
    Return the component if present, vector if not."""
    #slice(start,stop,step), None defaults to everything
    if True:
        #Centered difference: dv/dx_{i} = 1/(2 dx)(v_{i+1} - v_{i-1})
        SL = slice(None,-2,None)  #V_{i-1}
        SR = slice(2,None,None)   #V_{i+1}
        act=slice(1,-1,None)      #active region
        width = 2
    dx = nar(dx_in)*width
    output = [np.zeros_like( V[0]/dx_in[0]),np.zeros_like( V[0]/dx_in[0] ),np.zeros_like( V[0]/dx_in[0] ) ]

 
    #didj == dV_i/dx_j
    if component in [0,None]:
        d2d1 = 1./dx[1]*(V[2][act,SR,act] - V[2][act,SL,act])
        d1d2 = 1./dx[2]*(V[1][act,act,SR] - V[1][act,act,SL])
        output[0][act,act,act] = d2d1 - d1d2 
    if component in [1,None]:
        d0d2 = 1./dx[2]*(V[0][act,act,SR] - V[0][act,act,SL])
        d2d0 = 1./dx[0]*(V[2][SR,act,act] - V[2][SL,act,act])
        output[1][act,act,act] = d0d2 - d2d0 

    if component in [2,None]:
        d1d0 = 1./dx[0]*(V[1][SR,act,act] - V[1][SL,act,act])
        d0d1 = 1./dx[1]*(V[0][act,SR,act] - V[0][act,SL,act])
        output[2][act,act,act] = d1d0 - d0d1
    if component is not None:
        output = output[component]
    return output

def Cross(v1,v2):
    """Returns the cross product of v1 and v2.   More intuitive interface than numpy.cross"""
    return [v1[1]*v2[2] - v1[2]*v2[1],v1[2]*v2[0] - v1[0]*v2[2],v1[0]*v2[1] - v1[1]*v2[0]]

def _DivArray(array, dx, dy=None, dz=None,hydromethod=2): 
    # divergence of a vector field.
    if hydromethod == 2: 
        sl_left = slice(None,-2,None) 
        sl_right = slice(1,-1,None) 
        div_fac = 1.0 
    else: 
        sl_left = slice(None,-2,None) 
        sl_right = slice(2,None,None) 
        div_fac = 2.0 
    f  = array[0][sl_right,1:-1,1:-1]/dx 
    f -= array[0][sl_left ,1:-1,1:-1]/dx 
    if dy == None:
        dy = dx
    f += array[1][1:-1,sl_right,1:-1]/dy 
    f -= array[1][1:-1,sl_left ,1:-1]/dy 
    if dz == None:
        dz = dx
    f += array[2][1:-1,1:-1,sl_right]/dz 
    f -= array[2][1:-1,1:-1,sl_left ]/dz 
    new_field = np.zeros(array[0].shape, dtype='float64') 
    new_field[1:-1,1:-1,1:-1] = f 
    return new_field 

