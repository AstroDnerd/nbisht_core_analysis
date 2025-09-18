
from starter2 import *
import track_loader as TL
etrack_list=['m0230', 'm0231', 'm0232', 'm0233', 'm0234', 'm0235', 'm0236', 'm0237', 'm0238', 'm0239', 'm0240', 'm0241', 'm0242',\
          'm0243', 'm0244', 'm0245', 'm0246', 'm0247', 'm0250', 'm0260', 'm0270', 'm0280', 'm02100', 'm02110', 'u502', 'nb101']


TL.load_tracks(etrack_list)
def frame_from_et(name):
    if name.startswith('m'):
        return int(name[3:])
    else:
        return 118

class ucore_obj():
    def __init__(self,etname=None,particles=None,core_id=None, uid=None):
        self.uid=uid
        self.set_by_et={}
        self.core_id_by_et={}
        self.parents_by_et={}
        self.et_list=[]
        self.update(etname,this_set,core_id,[-1])
    def update(self,etname,this_set,core_id,parents):
        self.et_list.append(etname)
        self.set_by_et[etname] = this_set
        self.core_id_by_et[etname] = core_id
        self.parents_by_et[etname]=parents
        self.last_set = self.set_by_et[etname]


if 'ucore_list' not in dir():
    core_frame_ucore={}
    ucore_list=[]
    print('build ucore list')
    for etname in etrack_list:
        frame = frame_from_et(etname)
        this_looper=TL.tracks[etname]
        core_list = np.unique(this_looper.tr.core_ids)
        for core_id in core_list:
            print('=== core_id',core_id)
            if core_id not in core_frame_ucore:
                core_frame_ucore[core_id]={}
            particle_ids = np.unique(this_looper.tr.c([core_id],'particle_index'))
            this_set = set(particle_ids)
            keep_ucore = True
            noverlap=np.zeros(len(ucore_list))
            for nu,uc in enumerate(ucore_list):
                intersect = uc.last_set.intersection(this_set)
                if len(intersect) > 0:
                    noverlap[nu]=len(intersect)
            if noverlap.sum()>0:
                parents = np.arange(len(ucore_list))[nar(noverlap)>0]
                if (noverlap>0).sum()>1:
                    print(parents)
                ucore_id = np.argmax(noverlap)
                ouc = ucore_list[ucore_id]
                ouc.update(etname,this_set,core_id,parents)
                keep_ucore=False
            else:
                ucore_id = len(ucore_list)
                new_ucore=ucore_obj(etname,this_set,core_id,uid=ucore_id )
            core_frame_ucore[core_id][frame]=ucore_id

            if keep_ucore:
                ucore_list.append(new_ucore)
