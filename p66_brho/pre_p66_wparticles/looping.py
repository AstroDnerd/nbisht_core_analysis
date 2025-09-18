# %run -i loop_tools.py 
# For the command line after running loop_tools.py

from loop_tools import *

#ds = yt.load('/scratch1/dkl14b/add_gravity/GravPotential/DD0125/data0125')
ds = yt.load('/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/DD0125/data0125')
cores=get_leaf_indices(ds,h5_name='u05_0125_peaklist.h5')

coros = na.zeros(338,dtype=int)

for i in range(338):
    coros[i] = len(cores.get(i,0))

counter = 0
number = 0
fifteen = na.empty([0],dtype=int)

for i in range(338):
    if coros[i] > number:
         counter+=1 
         fifteen = na.append(fifteen,i)   
#        #print(i, end=' ')

#for i in range(338):
#    text = open("looping.text", "a")
#    text.write("\n Core#:%d, particles:%d"%(i,coros[i]))
#    text.close

#y = ",".join(map(str,fifteen))
#x = [y]
#text = open("looping.text", "a")
#text.write("Number of cores > %d particles: %d, Core numbers:\n %s" %(number,counter,x))
##text.write("\n %s" %coros)
#text.close




