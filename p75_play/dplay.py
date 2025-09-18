
from starter2 import *
import track_loader as TL
etrack_list=['m0230', 'm0231', 'm0232', 'm0233', 'm0234', 'm0235', 'm0236', 'm0237', 'm0238', 'm0239', 'm0240', 'm0241', 'm0242',\
          'm0243', 'm0244', 'm0245', 'm0246', 'm0247', 'm0250', 'm0260']#, 'm0270', 'm0280', 'm02100', 'm02110']  
#etrack_list = ['m0247','m0250','m0260']
#etrack_list = etrack_list[14:]#:17]
#etrack_list = etrack_list[:1]

def frame_from_et(name):
    return int(name[3:])

TL.load_tracks(etrack_list)

class ucore_obj():
    def __init__(self,etname=None,particles=None,core_id=None, uid=None):
        self.uid=uid
        self.set_by_et={}
        self.core_id_by_et={}
        self.parents_by_et={}
        self.et_list=[]
        self.update(etname,this_set,core_id,['start'])
    def update(self,etname,this_set,core_id,parents):
        self.et_list.append(etname)
        self.set_by_et[etname] = this_set
        self.core_id_by_et[etname] = core_id
        self.parents_by_et[etname]=parents
        self.last_set = self.set_by_et[etname]
    def get_ycoord(self):
        return self.uid


if 'ucore_list' not in dir() or False:
    core_frame_ucore={}
    ucore_list=[]
    for etname in etrack_list:
        frame = frame_from_et(etname)
        print(etname)
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
                print('new')
                ucore_id = len(ucore_list)
                new_ucore=ucore_obj(etname,this_set,core_id,uid=ucore_id )
            print('ucore',ucore_id)
            core_frame_ucore[core_id][frame]=ucore_id

            if keep_ucore:
                ucore_list.append(new_ucore)

if 0:
    ucore  = ucore_list[2]
    particles = list(ucore.last_set)
    this_looper = TL.tracks['m0246']
    for pid in particles:
        if pid in this_looper.tr.particle_ids:
            ind = np.where(this_looper.tr.particle_ids==pid)[0][0]
            core_id = this_looper.tr.core_ids[ind]
            print(core_id)

if 0:
    #plot the parent tree
    fig,axes=plt.subplots(1,1)
    ax0=axes

    for ucore in ucore_list[:20]:
        ycoord = ucore.get_ycoord()
        last_frame=-1
        for et in ucore.et_list:
            frame = frame_from_et(et)

            xcoord = frame

            ax0.scatter(xcoord, ycoord, c='r')
            if last_frame > 0:
                xcoord0 = last_frame
                if len(ucore.parents_by_et[et])>1:
                    for parent in ucore.parents_by_et[et]:
                        if frame != 60:
                            continue
                        print(frame,last_frame, ucore.uid)
                        ycoord0 = ucore_list[parent].get_ycoord()
                        coords=[[xcoord0,xcoord],[ycoord0,ycoord]]
                        ax0.plot(*coords)
                        print(coords)
            last_frame = frame

            #ax0.text(xcoord,ycoord,"%s"%ucore.uid)


    fig.savefig("%s/%s"%(plot_dir,'test'))


if 0:
    import tools.equal_probability_binner as epb
    ds0 = TL.tracks[etrack_list[0]].load(0)
    ad=ds0.all_data()
    rho0=ad['density']
    #hist0, cen0, width0 = epb.equal_prob(np.log10(rho0.v), 32)
    rhobins = np.geomspace(rho0.min(), rho0.max(),32)
    hist0, bins0 = np.histogram(rho0.v,bins=rhobins)
    bincen0 = 0.5*(bins0[1:]+bins0[:-1])
    for et in etrack_list:
        frame = frame_from_et(et)
        fig,axes=plt.subplots(1,2,figsize=(8,4))
        ax0=axes[0]; ax1=axes[1]
        this_looper = TL.tracks[et]
        core_list = np.unique(this_looper.tr.core_ids)
        stack_hist=0
        #for core_id in core_list:
        for ucore in ucore_list:
            if et not in ucore.core_id_by_et:
                continue
            core_id = ucore.core_id_by_et[et]
            ms = trackage.mini_scrubber(this_looper.tr,core_id,do_velocity=False)
            #x = ms.this_x.transpose()
            #y = ms.this_y.transpose()
            #ax.plot(x,y, c=[0.5]*4,linewidth=0.1)
            x = ms.mean_xc
            y = ms.mean_yc
            ax0.plot(x,y, c=[0.5]*4,linewidth=0.1)
            ax0.scatter(x[-1],y[-1],c='r',s=0.1)
            uid = core_frame_ucore[core_id][frame]
            ax0.text(x[-1],y[-1],"%d"%uid)

            this_rho = ms.density[:,0]
            hist, bins = np.histogram(this_rho,bins=rhobins)
            stack_hist = hist + stack_hist

            ax1.plot(bincen0, stack_hist)
        ax0.set(xlim=[-0.25,1.25],ylim=[-0.25,1.25])
        ax0.plot([0,1,1,0,0],[0,0,1,1,0],c='k')

        this_rho = this_looper.tr.track_dict['density'][:,0]
        hist, bins = np.histogram(this_rho,bins=rhobins)
        #ax1.plot(bincen0, hist0,c='k')
        ax1.plot(bincen0, hist,c='r')
        ax1.set(xscale='log',yscale='linear')
        #hist0, bins0 = np.histogram(rho0.v,bins=rhobins)
        #bincen0 = 0.5*(bins[1:]+bins[:-1])
        #ax1.bar(cen0,hist0,width=width0, color='k')
        #hist, cen, width = epb.equal_prob(np.log10(this_rho), 8)
        #ax1.bar(cen,hist,width=width, color='r')
        #ax1.set(yscale='log',xscale='linear')

        #for ucore in ucore_list:
        #    xcoord = 0.9
        #    ycoord = 0.98 - 0.05*ucore.uid
        #    ax0.text(xcoord,ycoord,"%d"%ucore.uid,transform=ax0.transAxes)



        fig.savefig('%s/%s'%(plot_dir,'test2_n%04d'%frame))

import convex_hull_tools as CHT
if 'ht' not in dir():
    ht={}
    for sim in etrack_list:
        ht[sim] = CHT.hull_tool(TL.loops[sim])
        ht[sim].make_hulls()

times=[]
for net,et in enumerate(etrack_list):
    times.append( TL.tracks[et].tr.times[-1])
if 1:
    ucore_volumes={}
    for ucore in ucore_list:
        ucore_volumes[ucore.uid]=np.zeros(len(etrack_list))
        for net,et in enumerate(etrack_list):
            if et in ucore.core_id_by_et:
                core_id = ucore.core_id_by_et[et]
            else:
                continue
            ind = np.where(ht[et].cores_used == core_id)[0][0]
            volume = ht[et].hull_volumes[ind]
            ucore_volumes[ucore.uid][net]=volume
    fig,ax=plt.subplots()
    cum=0
    for ucore in ucore_list:
        uid = ucore.uid
        cum = ucore_volumes[ucore.uid] + cum
        ax.plot(times/colors.tff,cum)
    ax.set(yscale='log')
    fig.savefig('%s/volume'%(plot_dir))




