from starter2 import *


skip_301 = [35,49]
skip_302 = [4, 89,17, 67, 109, 60, 58,59, 77,78, 90, 186, 188]
skip_302+= [100,102,302,303,310,313] #we should double check these.
skip_303 = [4,5,6,20] #these we have, but I think they're sus.
skip_303+= [10,122] #these are definitely problematic.

def verify_cores_u301(this_top):
    out = True
    if this_top.rhomax < 1100:
        out = False
    if this_top.peak_id in skip_301:
        out=False
    if this_top.leaf['particle_index'].size <10:
        out=False
    return out


def verify_cores_u302(this_top):
    out = True
    if this_top.rhomax < 1100:
        out = False
    if this_top.peak_id in skip_302:
        out=False
    if this_top.leaf['particle_index'].size <10:
        out=False
    return out

def verify_cores_u303(this_top):
    out = True
    if this_top.rhomax < 1100:
        out = False
    if this_top.peak_id in skip_303:
        out=False
    print('verify')
    if this_top.leaf['particle_index'].size <10:
        out=False
    return out

#things that are probably obsolete, but we'll keep until 
#we're sure we can reliably reproduce these results.

radius_u302 = {109:2e-2, 194:2e-2, 17:2e-2, 18:2e-2, 19:2e-2, 23:2e-2, 25:2e-2, 35:2e-2, 46:2e-2}
radius_u301 = {}
radius_u303 = {}

overlap_302 = {11: [14],
              #17: [18, 19],
              18: [19],
              20: [24, 25],
              21: [22],
              23: [25],
              26: [27, 29],
              27: [29],
              29: [30, 60],
              30: [60],
              33: [34, 35],
              34: [35, 46],
              35: [46],
              43: [45],
              45: [46],
              46: [52],
              50: [53, 56],
              52: [53],
              53: [55, 56, 57],
              55: [56, 57],
              56: [57],
              58: [59],
              60: [64],
              72: [73],
              77: [78],
              85: [86],
              87: [88, 95],
              90: [91, 92],
              91: [92],
              93: [104],
              100: [102],
              109: [110, 111],
              110: [111],
              186: [188],
              194: [195, 196],
              200: [201],
              208: [209],
              212: [213],
              273: [275, 277],
              275: [277],
              285: [286],
              302: [303],
              310: [313],
              349: [350],
              362: [363],
              366: [367, 368, 369, 370, 371, 372],
              367: [368, 369, 371, 372],
              368: [369, 371],
              369: [371, 372],
              371: [372],
              379: [380]}

overlap_303 = {4: [5, 6, 20],
             5: [6],
             10: [11, 12],
             11: [12],
             16: [17],
             17: [19],
             18: [19],
             19: [20],
             51: [52, 53],
             52: [53],
             93: [95, 96],
             95: [96],
             107: [108, 110],
             108: [110],
             111: [179],
             112: [113],
             115: [118],
             117: [118],
             121: [122, 123],
             122: [123],
             129: [139],
             132: [133, 135],
             133: [134, 135],
             134: [135],
             149: [150, 151, 152],
             150: [151],
             151: [152, 153, 155],
             153: [155],
             161: [167],
             170: [175, 179],
             172: [173],
             175: [179],
             177: [181],
             187: [188],
             194: [196, 198],
             196: [198],
             199: [201],
             204: [205],
             209: [211],
             235: [236],
             238: [239],
             270: [271]}
