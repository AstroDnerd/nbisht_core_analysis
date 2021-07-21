from starter2 import *

#coords = np.meshgrid( np.arange(0,10,0.5),np.arange(0,6,0.5) )
coords=np.mgrid[0:10:0.5,0:10:0.5]
x,y=np.mgrid[1:10:1,1:10:1]
coords = np.array([x,y])

myx,myy=coords[0].flatten(),coords[1].flatten()
coords_1=np.array(list(zip(myx,myy)))
F = np.array([[1.2,3],[0,1.0]])
if 0:
    this=coords
    X = np.einsum('ij,j...',F,this)
    X.shape
    print('all',coords)
    print('thi',this)
    print('pro',X)
if 1:
    b = np.array([ F@aa for aa in coords_1])
fig,axes=plt.subplots(2,1,figsize=(8,4),sharex=True)
axes[0].scatter(coords[0],coords[1])
bx = b[:,0]
by = b[:,1]
x1 = extents()
x1(x)
x1(bx)
y1 = extents()
y1(y)
y1(by)

axes[1].scatter(bx,by)
axbonk(axes[0],xlim=x1.minmax,ylim=y1.minmax)
axbonk(axes[1],xlim=x1.minmax,ylim=y1.minmax)

fig.savefig('plots_to_sort/cm_test.png')
plt.close('all')

if 'random1' not in dir():
    random1 =int((np.random.random()*by.size)//1)
    random2 =int((np.random.random()*by.size)//1)

dX = np.column_stack([b[random1],b[random2]])
dA = np.column_stack([coords_1[random1],coords_1[random2]])
print(dX)
print(dA)
Fp = np.linalg.solve(dA.T,dX.T).T
print('F prime')
print(Fp)
