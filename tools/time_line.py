"""
Take a figure, axis, and list of times/frames.
Draw a vertical line on the axis for each frame.
Save.
"""
from starter2 import *
import colors
def timelines(fig,ax,times,frames=None,name_template='PLOOT_%04d.png',**kwargs):

    for nt,time in enumerate(times):
        line = ax.axvline(time, **kwargs)
        outname = name_template%(frames[nt])
        fig.savefig(outname)
        print(outname)
        line.remove()



