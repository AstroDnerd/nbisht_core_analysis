from starter2 import *

# making 4 panel plots
fig,ax = plt.subplots(2,2, figsize=(6,6))
fig.subplots_adjust(wspace=0, hspace=0)

x = [0,1,2,3,4,5]
y = [0,3,4,6,8,10]
c = [3,4,5,6,7,8]
d = [6,8,9,7,9,2]

ax[0][0].scatter(x,y)
ax[0][1].scatter(y,x)
ax[1][0].scatter(c,d)
ax[1][1].scatter(y,x)

ax[0][0].xaxis.tick_top()
ax[0][1].xaxis.tick_top()
ax[0][1].xaxis.set_label_position('top')
ax[0][0].xaxis.set_label_position('top')
ax[1][1].yaxis.tick_right()
ax[0][1].yaxis.tick_right()
ax[1][1].yaxis.set_label_position('right')
ax[0][1].yaxis.set_label_position('right')

outname = 'test_four_way'
fig.savefig(outname)
