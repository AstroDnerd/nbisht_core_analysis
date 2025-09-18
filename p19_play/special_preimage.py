
from starter2 import *
import xtra_energy
import three_loopers_u500 as TL
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import pcolormesh_helper as pch
import colors
import movie_frames 
import other_scrubber
G = colors.G

import tsing
reload(tsing)
import radial
reload(radial)
import r_inflection
reload(r_inflection)
sim_list=['u501','u502','u503']
if 'tsing_tool' not in dir():
    tsing_tool={}
    for ns,sim in enumerate(sim_list):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()

def images(abbot, X, Y, subset=0):
    if subset==0:
        row_dict={}
        #row_dict['new_proj']=2
        row_dict['vt_cdf']=1
        row_dict['vr_cdf']=0
        row_dict['rho_cdf']=2

    Nplots=len(row_dict.keys())
    Ntimes=len(abbot.titles)
    extendor = [extents() for n in range(Nplots+1)]


    fig,axes=plt.subplots(Nplots,Ntimes, figsize=(12,12))
    for title,ax in zip(abbot.titles, axes[0]):
        ax.set_title(title)
    for core_id in abbot.cores_used:
        print('IMG',core_id)
        for nframe, frame in enumerate(abbot.frames[core_id]):
            row = row_dict['vr_cdf']
            ax = axes[row][nframe]
            ax.plot( X[core_id][frame]['vr_cdf_x_sph'], Y[core_id][frame]['vr_cdf_x_sph'], c='g')
            ax.plot( X[core_id][frame]['vr_cdf_x_par'], Y[core_id][frame]['vr_cdf_x_par'], c='purple')
            row = row_dict['vt_cdf']
            ax=axes[row][nframe]
            ax.plot( X[core_id][frame]['vt_cdf_x_sph'], Y[core_id][frame]['vt_cdf_x_sph'], c='g')
            ax.plot( X[core_id][frame]['vt_cdf_x_par'], Y[core_id][frame]['vt_cdf_x_par'], c='purple')
            row = row_dict['rho_cdf']
            ax=axes[row][nframe]
            ax.plot( X[core_id][frame]['rho_cdf_x_sph'], Y[core_id][frame]['rho_cdf_x_sph'], c='g')
            ax.plot( X[core_id][frame]['rho_cdf_x_par'], Y[core_id][frame]['rho_cdf_x_par'], c='purple')

    if 'rho_cdf' in row_dict:
        row=row_dict['rho_cdf']
        for ax in axes[row]:
            ax.set(xscale='log')

    outname = 'plots_to_sort/supper.png'
    fig.tight_layout()
    fig.savefig(outname)
def catering(abbot, suffix1='',subset=0):
    if subset==0:
        row_dict={}
        #row_dict['new_proj']=2
        row_dict['vt_cdf']=1
        row_dict['vr_cdf']=0
        row_dict['rho_cdf']=2

    Nplots=len(row_dict.keys())
    Ntimes=len(abbot.titles)
    extendor = [extents() for n in range(Nplots+1)]
    X={}
    Y={}


    fig,axes=plt.subplots(Nplots,Ntimes, figsize=(12,12))
    for title,ax in zip(abbot.titles, axes[0]):
        ax.set_title(title)
    for core_id in abbot.cores_used:
        X[core_id]={}
        Y[core_id]={}

        for iframe,frame in enumerate(abbot.frames[core_id]):
            X[core_id][frame]={}
            Y[core_id][frame]={}
            print('frame',frame)
            ds = abbot.ds[core_id][frame]
            reg = abbot.cubes[core_id][frame]
            sph = abbot.spheres[core_id][frame]
            R_KEEP = r_inflection.RKEEP(sph)
            ms = abbot.ms[core_id]
            nf = abbot.nf[core_id][frame]
            dv = sph[YT_cell_volume]
            RR = sph['radius']
            DD = sph[YT_density]
            EG = sph[YT_grav_energy_2]
            B2 = sph[YT_magnetic_field_strength]
            ORDER = np.argsort( RR)

            if 'vr_cdf' in row_dict:
                if 1:
                    if 0:
                        #probably the right thing to do
                        vel = []
                        for axis in 'xyz':
                            vel.append( sph['velocity_%s'%axis][ORDER][:10].mean())
                    if 1:
                        #the consistent thing to do
                        vel = ds.arr(ms.vcentral[:,nf], 'code_velocity')
                    scrub = other_scrubber.scrubber(sph, reference_velocity = vel)
                    scrub.compute_ke_rel()
                EK = scrub.ke_rel
                vt = scrub.vt_rel
                vr = scrub.vr_rel
                vr_particles = ms.vr_rel[:,nf]
                vt_particles = ms.vt2_rel[:,nf]**0.5
                rho_particles=ms.density[:,nf]

                row=row_dict['vr_cdf']
                ax=axes[row][iframe]

                cdf_vr = vr+0
                cdf_vr.sort()
                y_vr = np.arange(cdf_vr.size)/cdf_vr.size
                ax.plot(cdf_vr,y_vr, c='g')

                cdf_vp = vr_particles+0
                cdf_vp.sort()
                y_vp = np.arange(cdf_vp.size)/cdf_vp.size
                ax.plot(cdf_vp,y_vp, c='purple')
                extendor[row](cdf_vp)
                extendor[row](cdf_vr)
                X[core_id][frame]['vr_cdf_x_sph']=cdf_vr
                Y[core_id][frame]['vr_cdf_x_sph']=y_vr
                X[core_id][frame]['vr_cdf_x_par']=cdf_vp
                Y[core_id][frame]['vr_cdf_x_par']=y_vp

                row = row_dict['vt_cdf']
                ax=axes[row][iframe]
                cdf_vr = vt+0
                cdf_vr.sort()
                y_vr = np.arange(cdf_vr.size)/cdf_vr.size
                ax.plot(cdf_vr,y_vr, c='g')

                cdf_vp = vt_particles+0
                cdf_vp.sort()
                y_vp = np.arange(cdf_vp.size)/cdf_vp.size
                ax.plot(cdf_vp,y_vp, c='purple')
                extendor[row](cdf_vp)
                extendor[row](cdf_vr)
                X[core_id][frame]['vt_cdf_x_sph']=cdf_vr
                Y[core_id][frame]['vt_cdf_x_sph']=y_vr
                X[core_id][frame]['vt_cdf_x_par']=cdf_vp
                Y[core_id][frame]['vt_cdf_x_par']=y_vp


            if 'rho_cdf' in row_dict:
                row=row_dict['rho_cdf']
                ax=axes[row][iframe]
                rho_sphere=DD+0
                rho_sphere.sort()
                y_sphere = np.arange(rho_sphere.size)/rho_sphere.size
                rho_p = rho_particles+0
                rho_p.sort()
                y_p = np.arange(rho_p.size)/rho_p.size
                ax.plot(rho_sphere,y_sphere, c='g')
                extendor[row](rho_sphere)
                extendor[row](rho_p)
                ax.plot(rho_p,y_p, c='purple')
                X[core_id][frame]['rho_cdf_x_sph']=rho_sphere
                Y[core_id][frame]['rho_cdf_x_sph']=y_sphere
                X[core_id][frame]['rho_cdf_x_par']=rho_p
                Y[core_id][frame]['rho_cdf_x_par']=y_p


        if 'vr_cdf' in row_dict:
            row = row_dict['vr_cdf']
            for ax in axes[row]:
                xlim = np.abs(extendor[row].minmax).max()
                ax.set(xlim=[-xlim,xlim])
                ax.axhline(0.5)
                ax.axvline(0)
        if 'rho_cdf' in row_dict:
            row=row_dict['rho_cdf']
            for ax in axes[row]:
                ax.set(xscale='log',yscale='linear')
    outname = 'plots_to_sort/buffet_c%04d'%core_id

    fig.tight_layout()
    fig.savefig(outname)
    return X,Y




def delicatessen(abbot, suffix1='', subset=0):

    if subset==0:
        row_dict={'proj':0,'phi_phase':1,'rinf':-1}
        #row_dict['new_proj']=2
        row_dict['vr_cdf']=2
        row_dict['rho_cdf']=4

    Nplots=len(row_dict.keys())
    Ntimes=len(abbot.titles)
    extendor = [extents() for n in range(Nplots+1)]


    for core_id in abbot.cores_used:
        fig,axes=plt.subplots(Nplots,Ntimes, figsize=(12,12))
        for title,ax in zip(abbot.titles, axes[0]):
            ax.set_title(title)

        for iframe,frame in enumerate(abbot.frames[core_id]):
            print('frame',frame)
            ds = abbot.ds[core_id][frame]
            reg = abbot.cubes[core_id][frame]
            sph = abbot.spheres[core_id][frame]

            if 'proj' in row_dict:
                radius=sph.radius.v
                left = reg.left_edge.v
                right=reg.right_edge.v
                #center=reg.center.v
                center = 0.5*(left+right)

                width = (right-left)[0]
                proj_axis=0
                proj=yt.ProjectionPlot(ds,'xyz'[proj_axis],'density',center=center,data_source=reg)
                proj.set_width(width)
                image=proj.data_source.to_frb(width,128)['density']
                norm = mpl.colors.LogNorm(vmin=image[image>0].min(),vmax=image.max())
                row = row_dict['proj']
                ax = axes[row][iframe]

                L = list(left+0)
                R = list(right+0)
                C=list(center+0)
                L.pop(proj_axis)
                R.pop(proj_axis)
                C.pop(proj_axis)
                ext=[L[0],R[0],L[1],R[1]]

                ccc=C
                #pdb.set_trace()
                ax.imshow(image,norm=norm,extent=ext)
                circle= mpl.patches.Circle(ccc,radius=radius, ec='red',fc='None')
                ax.add_artist(circle)

            if 'phi_phase' in row_dict:
                row = row_dict['phi_phase']
                ax = axes[row][iframe]
                GE = np.abs(sph[YT_grav_energy_2])
                dv = np.abs(sph[YT_cell_volume])
                RR = sph[YT_radius]
                gbins = np.geomspace( GE[GE>0].min(), GE.max(),65)
                rbins = np.geomspace( RR [RR >0].min(), RR .max(),67)
                r_cen = 0.5*(rbins[1:]+rbins[:-1])
                hist, xb, yb = np.histogram2d( RR , GE, bins=[rbins,gbins],weights=dv)
                pch.helper(hist,xb,yb,ax=ax,transpose=False)
                ax.set(xscale='log',yscale='log')

            if 'rinf' in row_dict:
                #h2 is the histogram.
                #we'll remove any stragglers.
                h2 = hist+0
                shifter = np.zeros(nar(h2.shape)+2)
                cuml = np.zeros(h2.shape)
                c_center = slice(1,-1)
                #hb is just "is h2 nonzero"
                #we'll slide this around to look for neighbors
                hb = (h2>0)
                shifter[c_center,c_center]=hb
                nx,ny=shifter.shape
                for i in [0,1,2]:
                    for j in [0,1,2]:
                        if i==1 and j==1:
                            continue
                        s1 = slice(i,nx-2+i)
                        s2 = slice(j,ny-2+j)
                        cuml += shifter[s1,s2]
                #kill points that don't have neighbors.
                h2 *= (cuml >0)

                #Compute the upper bound of the histogram
                #smooth it
                #look for the point where the slope goes up.
                #but to avoid wiggles, it has to first come down a lot.

                #the upper bound of the distribution
                #compute upper_envelope
                y = np.arange(h2.shape[1])
                y2d = np.stack([y]*h2.shape[0])
                argmax = np.argmax(y2d*(h2>0),axis=1)
                upper_envelope = gbins[argmax]
                keepers=upper_envelope>1

                #smooth for wiggles.
                UE = gaussian_filter(upper_envelope[keepers],1)

                #
                # Find the inflection point where the slope comes up again.
                #
                #the slope
                DUE = UE[1:]-UE[:-1]
                #the max up to this radius
                cummax=np.maximum.accumulate( UE)
                #where the slope is positive
                ok = DUE > 0
                #and not too close to the center (for noise)
                ok = np.logical_and(ok , r_cen[keepers][1:]>1e-3)
                #and it has to come down below half its maximum
                ok = np.logical_and(ok, UE[1:]<0.5*cummax[1:])

                index = np.argmax(r_cen[keepers])
                if ok.any():
                    #find out the radius where the inflection happens.
                    index = np.where(ok)[0][0]
                R_KEEP = r_cen[keepers][index]
                print(R_KEEP)
                if iframe in [0,1]:
                    R_KEEP = radius


                circle= mpl.patches.Circle(ccc,radius=R_KEEP, ec='red',fc='None')
                ax = axes[ row_dict['proj']][iframe]
                ax.add_artist(circle)

                ax = axes[ row_dict['phi_phase']][iframe]

                ax.axvline(R_KEEP,c='r')
                R = abbot.ds[core_id][frame].arr(R_KEEP,'code_length')
                new_sphere = abbot.ds[core_id][frame].sphere(center,R)

            if 'new_proj' in row_dict:
                proj=yt.ProjectionPlot(ds,'xyz'[proj_axis],'density',center=center,data_source=new_sphere)
                proj.set_width(R_KEEP)
                image=proj.data_source.to_frb(width,128)['density']
                norm = mpl.colors.LogNorm(vmin=image[image>0].min(),vmax=image.max())
                row = row_dict['new_proj']
                ax = axes[row][iframe]

                L = list(left+0)
                R = list(right+0)
                C=list(center+0)
                L.pop(proj_axis)
                R.pop(proj_axis)
                C.pop(proj_axis)
                ext=[L[0],R[0],L[1],R[1]]

                ccc=C
                #pdb.set_trace()
                ax.imshow(image,norm=norm,extent=ext)

            if 'vr_cdf' in row_dict:
                print('word')
                ms = abbot.ms[core_id]
                nf = abbot.nf[core_id][frame]
                dv = sph[YT_cell_volume]
                RR = sph['radius']
                DD = sph[YT_density]
                EG = sph[YT_grav_energy_2]
                B2 = sph[YT_magnetic_field_strength]

                ORDER = np.argsort( RR)
                if 1:
                    if 0:
                        #probably the right thing to do
                        vel = []
                        for axis in 'xyz':
                            vel.append( sph['velocity_%s'%axis][ORDER][:10].mean())
                    if 1:
                        #the consistent thing to do
                        vel = ds.arr(ms.vcentral[:,nf], 'code_velocity')
                    scrub = other_scrubber.scrubber(sph, reference_velocity = vel)
                    scrub.compute_ke_rel()
                EK = scrub.ke_rel
                vt = scrub.vt_rel
                vr = scrub.vr_rel
                vr_particles = ms.vr_rel[:,nf]
                vt_particles = ms.vt2_rel[:,nf]**0.5
                rho_particles=ms.density[:,nf]

                row=row_dict['vr_cdf']
                ax=axes[row][iframe]

                cdf_vr = vr+0
                cdf_vr.sort()
                y_vr = np.arange(cdf_vr.size)/cdf_vr.size
                ax.plot(cdf_vr,y_vr, c='g')

                cdf_vp = vr_particles+0
                cdf_vp.sort()
                y_vp = np.arange(cdf_vp.size)/cdf_vp.size
                ax.plot(cdf_vp,y_vp, c='purple')
                extendor[row](cdf_vp)
                extendor[row](cdf_vr)

                ax=axes[row+1][iframe]
                cdf_vr = vt+0
                cdf_vr.sort()
                y_vr = np.arange(cdf_vr.size)/cdf_vr.size
                ax.plot(cdf_vr,y_vr, c='g')

                cdf_vp = vt_particles+0
                cdf_vp.sort()
                y_vp = np.arange(cdf_vp.size)/cdf_vp.size
                ax.plot(cdf_vp,y_vp, c='purple')
                extendor[row](cdf_vp)
                extendor[row](cdf_vr)

            if 'rho_cdf' in row_dict:
                row=row_dict['rho_cdf']
                ax=axes[row][iframe]
                rho_sphere=DD+0
                rho_sphere.sort()
                y_sphere = np.arange(rho_sphere.size)/rho_sphere.size
                rho_p = rho_particles+0
                rho_p.sort()
                y_p = np.arange(rho_p.size)/rho_p.size
                ax.plot(rho_sphere,y_sphere, c='g')
                extendor[row](rho_sphere)
                extendor[row](rho_p)
                ax.plot(rho_p,y_p, c='purple')


        if 'vr_cdf' in row_dict:
            row = row_dict['vr_cdf']
            for ax in axes[row]:
                xlim = np.abs(extendor[row].minmax).max()
                ax.set(xlim=[-xlim,xlim])
                ax.axhline(0.5)
                ax.axvline(0)
        if 'rho_cdf' in row_dict:
            row=row_dict['rho_cdf']
            for ax in axes[row]:
                ax.set(xscale='log',yscale='linear')
        outname = 'plots_to_sort/samich_c%04d'%core_id

        fig.tight_layout()
        fig.savefig(outname)








    



class butcher():
    def __init__(self, this_looper):
        self.this_looper=this_looper

        self.spheres={}
        self.cubes={}
        self.ms={}
        self.cores_used=[]
        self.frames={}
        self.nf={}
        self.ds={}

    def run(self, core_list=None, tsing=None, square_factor=1, timescale=0):

        this_looper=self.this_looper
        thtr=this_looper.tr
        self.square_factor=square_factor
        self.timescale=timescale
        def get_time_index(time):
            index=np.argmin( np.abs( thtr.times/colors.tff-time))
            return index

        for core_id in core_list:
            frame_mask = np.zeros_like(thtr.times, dtype='bool')

            if self.timescale==0:
                frame_mask[0]=True
                frame_mask[get_time_index(0.25*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(0.5*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(0.75*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(tsing.tend_core[core_id])]=True
                #half_collapse = 0.5*(tsing.tsing_core[core_id]+tsing.tend_core[core_id])
                #theta=0.5
                #half_collapse = theta*tsing.tsing_core[core_id]+(1-theta)*tsing.tend_core[core_id]
                #frame_mask[get_time_index(half_collapse)]=True
                self.titles=[ r'$t=0$', r'$t=0.25 t_{\rm{sing}}$', r'$t=0.5 t_{\rm{sing}}$', r'$t=0.75 t_{\rm{sing}}$',\
                    r'$t=t_{\rm{sing}}$', r'$t=t_{\rm{sung}}$']
            if self.timescale==1:
                self.titles=[]
                for theta in np.linspace(0,1,Ntimes):
                    t = (1-theta)*tsing.tsing_core[core_id]+theta*tsing.tend_core[core_id]
                    self.titles.append(r'%0.2f'%theta)
                    index=get_time_index(t)
                    frame_mask[index]=True
            if self.timescale==2:
                frame_mask[0]=True
                #frame_mask[get_time_index(0.25*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(0.5*tsing.tsing_core[core_id])]=True
                #frame_mask[get_time_index(0.75*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(tsing.tend_core[core_id])]=True
                #half_collapse = 0.5*(tsing.tsing_core[core_id]+tsing.tend_core[core_id])
                #theta=0.5
                #half_collapse = theta*tsing.tsing_core[core_id]+(1-theta)*tsing.tend_core[core_id]
                #frame_mask[get_time_index(half_collapse)]=True
                self.titles=[ r'$t=0$', r'$t=0.5 t_{\rm{sing}}$', r'$t=t_{\rm{sing}}$', r'$t=t_{\rm{sung}}$']
            if self.timescale==3:
                frame_mask[get_time_index(tsing.tend_core[core_id])]=True
                self.titles=[  r'$t=t_{\rm{sung}}$']
            if sum(frame_mask) == 0:
                pdb.set_trace()
            
            self.frames[core_id]=thtr.frames[frame_mask]
            self.spheres[core_id]={}
            self.cubes[core_id]={}
            self.ds[core_id]={}
            self.nf[core_id]={}
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            ms.get_central_at_once(core_id)
            ms.compute_ge(core_id)
            ms.compute_ke_rel(core_id)
            self.ms[core_id]=ms
            self.cores_used.append(core_id)

            for nframe,frame in enumerate(self.frames[core_id]):
                ds = this_looper.load(frame)
                self.ds[core_id][frame]=ds
                xtra_energy.add_energies(ds)
                xtra_energy.add_gdotgradrho(ds)
                nf = np.where( this_looper.tr.frames == frame)[0][0]
                self.nf[core_id][frame]=nf

                center = nar([ms.mean_xc[nf], ms.mean_yc[nf],ms.mean_zc[nf]])
                msR = ms.rc
                msR[ msR<1/2048]=1/2048

                MaxRadius=msR[:,nf].max()
                Radius = max([8.0/128, MaxRadius])
                rsph = ds.arr(Radius,'code_length')
                self.spheres[core_id][frame] = ds.sphere(center,rsph)

                left = center - Radius*self.square_factor
                right = center + Radius * self.square_factor
                center = 0.5*(left+right)

                self.cubes[core_id][frame]=ds.region(center,left,right)








sim_list=['u501','u502','u503']
sim_list=['u501']
if 'abattoir' not in dir():
    abattoir={}
for sim in sim_list:
    if sim not in abattoir:
        all_cores=np.unique( TL.loops[sim].tr.core_ids)
        core_list=list(all_cores)
        core_list=None#[323]
        core_list=[323]
        core_list=[25]
        #core_list=[74]
        core_list = TL.loops[sim].core_by_mode['Alone']
        #core_list=core_list[:2]
        this=butcher(TL.loops[sim])
        timescale = 2 #0= 0-tsing, 1=tsing-tsing 2=4 panel
        this.run(core_list=core_list,tsing=tsing_tool[sim], timescale=timescale)#, r_inflection=anne.inflection[sim])
        abattoir[sim]=this
if 0:
    for sim in sim_list:
        #phasor(mp[sim])
        delicatessen(abattoir[sim])
if 1:
    for sim in sim_list:
        #phasor(mp[sim])
        if 'ThisX' not in dir():
            ThisX, ThisY=catering(abattoir[sim])
        images(abattoir[sim], ThisX, ThisY)
