from starter1 import *

def add_sublabel(ax,input_tool):
    """Apply the sublable.  
    Might be hung onto a tool, might be free standing."""
    if hasattr(input_tool,'sublabel'):
        if input_tool.sublabel is not None:
            input_tool.sublabel(ax)
    elif type(input_tool) == labs:
        input_tool(ax)

    

class labs():
    def __init__(self,xpos,ypos,label,**kwargs):
        self.xpos=xpos
        self.ypos=ypos
        self.label=label
        self.kwargs=kwargs
    def __call__(self,ax):
        ax.text( self.xpos, self.ypos, self.label, **self.kwargs)
