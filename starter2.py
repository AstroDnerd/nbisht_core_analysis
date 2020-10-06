
#
# ALl the other things that aren't standard python packages.
#

from starter1 import *
#so we can import things from sub directories
path_list= ["./tools_data", "./tools_pdf", "./testing",\
            "./tools", "./tools_tracks", "./trash"]
for directory in path_list:
    if directory not in sys.path:
        sys.path += [directory]

import yt
import looper
reload(looper)
import trackage
reload(trackage)
#from core_dump import *

import loop_tools
reload(loop_tools)

import looper
reload(looper)
import trackage
reload(trackage)
import tracks_read_write as trw
reload(trw)
from davetools import *
import loop_apps
reload(loop_apps)



