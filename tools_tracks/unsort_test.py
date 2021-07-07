

N=10
a = nar([int( (100*np.random.random())//1) for i in range(10)])
asort = np.argsort(a)
a1=a[asort]
rev = np.argsort(asort)
a1[rev]
print(a)
print(a1)
print(a1[rev])
