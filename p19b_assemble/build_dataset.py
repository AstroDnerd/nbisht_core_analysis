from dtools.starter1 import *
import yt



def write_enzo_set(arr_in, directory,field=None, do_swap=True, dims=None, extend=None):

    if do_swap:
        arr_in = arr_in.swapaxes(0,2)
    if extend is not None:
        m_extend = 2-extend
        #m_extend = extend
        shape = nar(arr_in.shape)
        shape[m_extend]+=1
        all_slice = [slice(None),slice(None),slice(None)]
        cube = tuple([slice(0,arr_in.shape[0]),
                slice(0,arr_in.shape[1]),
                slice(0,arr_in.shape[2])])
        face1 = copy.copy(all_slice)
        face2 = copy.copy(all_slice)
        face1[m_extend]=0
        face2[m_extend]=-1
        face1=tuple(face1)
        face2=tuple(face2)
        arr = np.zeros(shape)
        arr[cube]=arr_in
        arr[face2]=arr_in[face1]
    else:
        arr=arr_in
    print(field, arr.shape)
    




    print('Writing',field)
    fname = "%s/%s"%(directory,field)
    fptr=h5py.File(fname,'w')
    fptr[field]=np.ascontiguousarray(arr)
    fptr[field].attrs['Component_Rank']=1
    fptr[field].attrs['Component_Size']=arr.size #this value is actually irrelevant.
    fptr[field].attrs['Rank']=len(arr.shape)
    if dims is None:
        dims = arr.shape
    fptr[field].attrs['Dimensions']=dims[::-1]

    fptr.close()

