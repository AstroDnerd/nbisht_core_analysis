from starter1 import *

def unique(L):
    """A terrible n^2 search for uniqueness in a list.
    This list might contain tuples, because yt."""
    output=[]
    for value in L:
        got_it=False
        for o in output:
            if value == o:
                got_it=True
                break
        if not got_it:
            output.append(value)
    return output
        
