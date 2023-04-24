from numpy import cross, eye, dot
from scipy.linalg import expm, norm

def M(axis, theta):
    return expm(cross(eye(3), axis/norm(axis)*theta))

v, axis, theta = nar([1,1,7]), nar([0,1,0]), 0.8
v=v/(v*v).sum()**0.5
M0 = M(axis, theta)
RV=dot(M0,v)
print(v2)
theta_w=np.arccos(dot(v,RV)/((v*v).sum()*(RV*RV).sum())**0.5)
print('one',(v*v).sum())
print('wtf',theta_w, np.arccos(np.cos(theta)))
# [ 2.74911638  4.77180932  1.91629719]

ahat = axis/(axis*axis).sum()
v1 = (ahat*v).sum()*ahat
v2 = v-v1
v2mag = np.sqrt((v2*v2).sum())
Dprime = v2mag*np.sin( theta/2)
vmag = np.sqrt((v*v).sum())
theta2 = 2*np.arcsin( Dprime/vmag)
#print(np.sin(theta_w/2), v2mag/vmag*np.sin(theta/2))
print("want %0.3f got %0.3f"%(theta_w, theta2))

RV2 = dot(M0,v2)
v3 = v1+RV2
#print(v3)
#print(RV)
#print(RV2)

Drot = v2-RV2
Dprime_2 = (Drot*Drot).sum()**0.5
print(Dprime_2/2, Dprime)

