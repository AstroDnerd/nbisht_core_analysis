from starter2 import *

import convex_hull_tools as CHT
import matplotlib.colors as colors

reload(CHT)
import hair_dryer
reload(hair_dryer)
import stay_close
import three_loopers_tenfour as TL4
sim_list=['u401','u402','u403']
#sim_list=['u402']
if 'ht' not in dir() :
    ht = {}
    for this_simname in sim_list:
        ht[this_simname] = CHT.hull_tool(TL4.loops[this_simname])
        ht[this_simname].make_hulls()
        ht[this_simname].make_overlaps()

if 'ct' not in dir():
    ct = {}
    for this_simname in sim_list:
        ct[this_simname] = stay_close.close_tool( TL4.loops[this_simname])
        ct[this_simname].make_distance()

import supersets
reload(supersets)
if 'st' not in dir():
    st={}
    for this_simname in sim_list:
        st[this_simname] = supersets.superset( TL4.loops[this_simname], ht[this_simname])
        st[this_simname].find()

if 'overlap_dict' not in dir():
    overlap_dict={}
    overlap_numb={}
    particles = {}
    ratio_matrix={}
    ratio_matrixb={}
    particle_matrix={}

    for ns,this_simname in enumerate(sim_list):
        htool = ht[this_simname]
        overlap_dict[this_simname] = np.zeros( [len(htool.cores_used)]*2) -1
        overlap_numb[this_simname] = np.zeros( [len(htool.cores_used)]*2) -1
        ratio_matrix[this_simname] = np.zeros( [len(htool.cores_used)]*2) 
        ratio_matrixb[this_simname] = []
        particle_matrix[this_simname] = np.zeros( [len(htool.cores_used)]*2)
        particles[this_simname]=[]
        for nc1,core_id_1 in enumerate(htool.cores_used):
            particles[this_simname].append(htool.nparticles[nc1])
        for nc1,core_id_1 in enumerate(htool.cores_used):
            for nc2,core_id_2 in enumerate(htool.cores_used):
                val = htool.overlaps[core_id_1][nc2]
                num = htool.overlap_number[core_id_1][nc2]
                overlap_dict[this_simname][nc1,nc2] = val
                overlap_numb[this_simname][nc1,nc2] = num

                a,b=htool.overlaps[core_id_1][nc2], htool.overlaps[core_id_2][nc1] 
                ratio = max([a,b])
                rat= sorted( [a,b])
                if rat[1] == 0:
                    ratio = max([a,b])
                else:
                    ratio=rat[0]/rat[1]
                ratio_matrix[this_simname][nc1,nc2]=ratio
                particle_matrix[this_simname][nc1,nc2] = np.sqrt( particles[this_simname][nc1]*particles[this_simname][nc2])



import means_etc
reload(means_etc)

if 1:
    #Probably the good one.
    for ns,this_simname in enumerate(sim_list):
        on = overlap_numb[this_simname]
        ont = on.transpose()
        of = overlap_dict[this_simname]
        oft = overlap_dict[this_simname].transpose()
        #how it works:
        #tool.overlaps[core_id] = list of foreign particles in my contour
        #o[core_id,:] = tool.overlaps[core_id]
        #o[core,:] = other particle in my core contour
        #o[:,core] = other contours with my particles
        both = np.stack([on, ont])
        minmin = both.min(axis=0)

        ofmin = np.stack([of,oft]).min(axis=0)
        of = ofmin


        mask = (on>0.0)*(ont>0.0)
        overlap = ofmin>0

        m1=mask*overlap
        NumberOfOverlap=(m1).sum(axis=1)
        Max = (of*mask).max(axis=1)
        NumberOfOverlap=(m1).sum(axis=1)
        Fraction = (of*mask).sum(axis=1)
        Fraction[NumberOfOverlap>0] /= NumberOfOverlap[ NumberOfOverlap > 0]

        nmax = NumberOfOverlap.max()
        bins_n=np.linspace(-0.5,nmax+0.5)
        bins_f = np.linspace(0,1,20)

        plt.clf()
        figa, axa, axtop,axright = means_etc.three_way_bean()
        axa.scatter( NumberOfOverlap, Fraction,c='k')
        axa.scatter( NumberOfOverlap,Max ,c='r')
        axtop.hist( NumberOfOverlap, bins=bins_n, histtype='step',color='k')
        axright.hist( Fraction, bins=bins_f, histtype='step',color='k',orientation='horizontal')
        axright.hist( Max, bins=bins_f, histtype='step',color='r',orientation='horizontal')

        if 1:
            axbonk(axa,xlabel=r'$N_{\rm{overlap}}$', ylabel='Overlap Fraction')
            axa.set_xlim([-0.1,nmax+0.1])
            axbonk(axtop,xlabel='',ylabel=r'$N$')
            axbonk(axright,xlabel=r'$N$',ylabel='')
            axright.set_ylim( axa.get_ylim())
            axtop.set_xlim( axa.get_xlim())
            axtop.set_xticks([])
            axright.set_yticks([])

        figa.savefig('plots_to_sort/%s_overlaps.pdf'%this_simname)
if 0:
    fig,ax=plt.subplots(1,2)#, figsize=(12,4))
    for ns,this_simname in enumerate(sim_list):
        of = overlap_dict[this_simname]
        oft = of.transpose()
        o1 = of.flatten()
        o2 = oft.flatten()

        no = ~((o1==0)*(o2==0))
        no = no *( (o1 >=0)+(o2>=0))
        print((o1[no]==0).sum())
        b = np.column_stack([o1[no],o2[no]])
        c = np.min(b,axis=1)
        ax[0].hist(c,histtype='step',label=this_simname,bins=64)

        no = (( o1 > 0) + (o2 > 0) )*( ~( (o1 ==0)+(o2==0)))

    axbonk(ax[0],xlabel=r'$r1/r2$',ylabel='N',yscale='log')
    fig.savefig('plots_to_sort/hist.png')

if 0:
    fig,ax=plt.subplots(1,2)#, figsize=(12,4))
    for ns,this_simname in enumerate(sim_list):
        of = overlap_dict[this_simname]
        oft = of.transpose()
        o1 = of.flatten()
        o2 = oft.flatten()

        ok = ~((o1==0)*(o2==0))
        ok = ok *( (o1 >=0)+(o2>=0))
        print((o1[ok]==0).sum())
        b = np.column_stack([o1[ok],o2[ok]])
        c = np.min(b,axis=1)
        geo = np.geomspace(1e-2,1,16)
        bins = np.concatenate( [[0],geo])
        ax[0].hist(c,histtype='step',label=this_simname,bins=bins, density=True)

        xor = np.logical_xor( (o1==0), (o2==0))
        bins = np.geomspace( 1e-5,1,16)
        butins=o1[xor]
        butins=butins[butins>0]
        ax[1].hist(butins, histtype='step',label=this_simname,bins=bins, density=True)
        #ax[1].plot( sorted(butins))


    axbonk(ax[0],xlabel=r'$r1/r2$',ylabel='N',yscale='log',xscale='log')
    axbonk(ax[1],xlabel=r'Lopsides',ylabel='N',yscale='log',xscale='log')
    fig.savefig('plots_to_sort/hist.png')
        
        




if 0:

    fig,ax=plt.subplots(1,3, figsize=(12,4))
    for ns,this_simname in enumerate(sim_list):
        on = overlap_numb[this_simname]
        ont = on.transpose()
        of = overlap_dict[this_simname]
        #how it works:
        #tool.overlaps[core_id] = list of foreign particles in my contour
        #o[core_id,:] = tool.overlaps[core_id]
        #o[core,:] = other particle in my core contour
        #o[:,core] = other contours with my particles

        mask = (on>0.0)*(ont>0.0)

        part = particle_matrix[this_simname].flatten()
        ratt = ratio_matrix[this_simname].flatten()

        ax[ns].scatter(part,ratt)
        overlap = on>0
        m1=mask*overlap
        NumberOfOverlap=(m1).sum(axis=1)
        Fraction = (of*mask).sum(axis=1)
        Fraction[NumberOfOverlap>0] /= NumberOfOverlap[ NumberOfOverlap > 0]
    fig.savefig('plots_to_sort/part_rat.png')





if 0:
    fig,ax=plt.subplots(1,3)
    fig2,ax2=plt.subplots(1,3, figsize=(12,4))
    for ns,this_simname in enumerate(sim_list):
        on = overlap_numb[this_simname]
        ont = on.transpose()
        of = overlap_dict[this_simname]
        #how it works:
        #tool.overlaps[core_id] = list of foreign particles in my contour
        #o[core_id,:] = tool.overlaps[core_id]
        #o[core,:] = other particle in my core contour
        #o[:,core] = other contours with my particles

        mask = (on>0.0)*(ont>0.0)
        overlap = on>0
        this_in_that=overlap.sum(axis=0)
        that_in_this=overlap.sum(axis=1)
        nmax = max([this_in_that.max(),that_in_this.max()])
        bins=np.linspace(-0.5,nmax+0.5)

        ax[ns].hist(this_in_that, histtype='step', label='this in that',bins=bins)
        ax[ns].hist(that_in_this, histtype='step', label='that in this',bins=bins)

        m1=mask*overlap
        NumberOfOverlap=(m1).sum(axis=1)
        Fraction = (of*mask).sum(axis=1)
        Fraction[NumberOfOverlap>0] /= NumberOfOverlap[ NumberOfOverlap > 0]


        ax2[ns].scatter( NumberOfOverlap, Fraction,s=1)
        ax2[ns].set_title('No overlap %d'%(NumberOfOverlap==0).sum())

        if 0:
            bx_min=-0.5; bx_max=NumberOfOverlap.max()+1.5
            dx=1
            by_min=0; by_max=1.01; dy=0.01; Ny=int(1/dy)+1
            bins_x = np.arange(bx_min,bx_max,dx)
            bins_y = np.linspace(by_min,by_max,Ny)
            dy_actual = bins_y[1]-bins_y[0]
            h, xb,yb = np.histogram2d(NumberOfOverlap, Fraction,bins=(bins_x,bins_y))
            norm = colors.LogNorm(vmin=h[h>0].min(),vmax=h.max())
            cmap=copy.copy(mpl.cm.get_cmap("viridis"))
            #cmap.set_under('w')
            ix = np.floor(NumberOfOverlap).astype('int')
            iy = np.floor(Fraction/dy_actual).astype('int')
            vals = h[ix,iy]
            print(h[0,0])
            #print(vals)
            c=cmap(norm(vals))
            norms=c.sum(axis=1)
            #print( bins_x[ ix] - NumberOfOverlap)
            #pdb.set_trace()


            ax2[ns].scatter(NumberOfOverlap,Fraction,cmap=cmap,norm=norm,c=vals,s=1)
            #plot=ax2[ns].pcolormesh(the_x,the_y,h,norm=norm, cmap=cmap)

        if 0:
            bins_x = np.arange(-0.5,41.5)
            bins_y = np.linspace(0,1,10)
            h, xb,yb = np.histogram2d(NumberOfOverlap, Fraction,bins=(bins_x,bins_y))
            xc = 0.5*(xb[1:]+xb[:-1])
            yc = 0.5*(yb[1:]+yb[:-1])
            the_x = np.r_[ (yc.size)*[xc]].transpose()
            the_y = np.r_[ (xc.size)*[yc]]
            norm = colors.LogNorm()#vmin=1,vmax=h.max())
            cmap=copy.copy(mpl.cm.get_cmap("viridis"))
            cmap.set_under('w')
            plot=ax2[ns].pcolormesh(the_x,the_y,h,norm=norm, cmap=cmap)


    ax[0].legend(loc=0)
    fig.savefig('plots_to_sort/flail.png')
    fig2.savefig('plots_to_sort/flail2.png')
    plt.close(fig)


if 0:
    fig,ax=plt.subplots(1,3)
    for ns,this_simname in enumerate(sim_list):
        o = overlap_dict[this_simname]

        p1 = o.flatten()
        p2 = o.transpose().flatten()
        ok = (p1>-0.5)*(p2>-0.5) 
        ok[ (p1<0.01)+(p2<0.01)]=False
        p1b=p1[ok]
        p2b=p2[ok]
        #ax[ns].scatter( p1b,p2b)
        h, xb,yb = np.histogram2d(p1b,p2b,bins=10)
        xc = 0.5*(xb[1:]+xb[:-1])
        yc = 0.5*(yb[1:]+yb[:-1])
        the_x = np.r_[ (yc.size)*[xc]].transpose()
        the_y = np.r_[ (xc.size)*[yc]]
        norm = colors.LogNorm(vmin=1,vmax=h.max())
        cmap=copy.copy(mpl.cm.get_cmap("viridis"))
        cmap.set_under('w')
        plot=ax[ns].pcolormesh(the_x,the_y,h,norm=norm, cmap=cmap)

        fig3d=plt.figure()
        ax3d=fig3d.add_subplot(projection='3d')
        ax3d.plot_wireframe( the_x,the_y,h)
        fig3d.savefig('plots_to_sort/%s_overlaps.png'%this_simname)

    fig.savefig('plots_to_sort/test.png')
    plt.close(fig)

if 0:
    fig,ax=plt.subplots(1,3)
    for ns,this_simname in enumerate(sim_list):
        o = overlap_dict[this_simname]
        p1 = o.flatten()
        p2 = o.transpose().flatten()
        ok = (p1>-0.5)*(p2>-0.5) 
        #ok[ (p1<0.01)+(p2<0.01)]=False
        p1b=p1[ok]
        p2b=p2[ok]
        top = np.column_stack([p1b,p2b]).min(axis=1)
        bot = np.column_stack([p1b,p2b]).max(axis=1)
        #ax[ns].hist( top/bot)
        ax[ns].hist( top)
        ax[ns].set_yscale('log')
        #ax[ns].scatter( p1b,p2b)
    fig.savefig('plots_to_sort/overlap_ratio.png')

if 0:
    #every pair
    import means_etc
    for this_simname in sim_list:

        htool = ht[this_simname]
        ctool = ct[this_simname]
        stool = st[this_simname]

        hair = hair_dryer.hair_tool( htool.this_looper)
        overlap_matrix = np.zeros( [len(htool.cores_used)]*2) -1
        for nc1,core_id_1 in enumerate(htool.cores_used):
            for nc2,core_id_2 in enumerate(htool.cores_used):
                #if nc2 <= nc1:
                #    continue
                overlap_matrix[nc1,nc2] = htool.overlaps[core_id_1][nc2]
                continue


                if htool.overlaps[core_id_1][nc2] == 0 and htool.overlaps[core_id_2][nc1] > 0:
                    #core_1=core_id_1
                    #core_2=core_id_2
                    #self=htool
#
#                    hull_1 =  self.hulls[core_1]
#                    hull_2 =  self.hulls[core_2]
#                    vert_1 = self.points_3d[core_1][hull_1.vertices,:]
#                    vert_2 = self.points_3d[core_2][hull_2.vertices,:]
#                    points_1 = self.points_3d[core_1]
#                    points_2 = self.points_3d[core_2]
#
#                    in_1_2_3 = CHT.in_hull(points_1, points_2)
#                    fraction =  in_1_2_3.sum()/points_1.shape[0]


                    cdict={core_id_1:'r',core_id_2:'g'}
                    CHT.plot_watershed(htool, core_list=[core_id_1,core_id_2],frames=[0],accumulate=True, 
                                prefix='c%04d_c%04d'%(core_id_1,core_id_2), axis_to_plot=[-1],
                                color_dict=cdict,label_cores=[-1])
                    #hair.run( core_list=[core_id_1,core_id_2], newplots=False, colors=cdict, name = "c%04d_c%04d_"%(core_id_1,core_id_2))
            else:
                continue
            break




