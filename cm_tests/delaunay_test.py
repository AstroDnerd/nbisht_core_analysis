from starter2 import *
#https://stackoverflow.com/questions/50517812/find-all-simplices-a-point-is-a-part-of-in-scipy-spatial-delaunay-python

points = np.array([[0, 0], [0, 1.1], [1, 0], [1, 1]])
from scipy.spatial import Delaunay
tri = Delaunay(points)
def all_simplices(d,nvert):
    return [i for i,s in enumerate(d.simplices) if nvert in s]

print(all_simplices(tri,0))

