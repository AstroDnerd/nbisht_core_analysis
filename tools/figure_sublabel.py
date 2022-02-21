from starter1 import *

def add_sublabel(ax,looper):
    if hasattr(looper,'sublabel'):
        if looper.sublabel is not None:
            looper.sublabel(ax)
    

class labs():
    def __init__(self,xpos,ypos,label,**kwargs):
        self.xpos=xpos
        self.ypos=ypos
        self.label=label
        self.kwargs=kwargs
    def __call__(self,ax):
        ax.text( self.xpos, self.ypos, self.label, **self.kwargs)
