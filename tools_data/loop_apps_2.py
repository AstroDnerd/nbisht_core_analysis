
@looper.core_loop
def core_hist(looper,snapshot, field='density', axis_list=[0,1,2], color='r',force_log=None,linthresh=100):
    if snapshot.R_centroid is None:
        snapshot.get_all_properties()
    ds = snapshot.get_ds()
    for ax in axis_list:
        center = ds.arr(snapshot.R_centroid,'code_length')
        Rmax = snapshot.R_mag.max()
        scale_min = ds.arr(0.05,'code_length')
        scale = max([Rmax, scale_min])
        sph = ds.sphere(center,scale)
        prof = yt.create_profile(sph,field,'cell_volume')

from scipy.optimize import curve_fit
def powerlaw(r,rho0, r0, alpha):
    return alpha*np.log10(r/r0) + np.log10(rho0)
@looper.core_loop
def core_radial(looper,snapshot, fields=[],weight_field=None,scales=['log','log'] ):

    if "alpha_dict_part" not in looper.__dict__:
        looper.alpha_dict_part={}
        looper.alpha_dict_prof={}
    if looper.current_frame not in looper.alpha_dict_part:
        looper.alpha_dict_part[looper.current_frame]={}
        looper.alpha_dict_prof[looper.current_frame]={}

    plt.clf()
    if snapshot.R_centroid is None:
        snapshot.get_all_properties()
    ds = snapshot.get_ds()
    center = ds.arr(snapshot.R_centroid,'code_length')
    Rmax = snapshot.R_mag.max()
    scale_min = ds.arr(0.05,'code_length')
    scale = max([Rmax, scale_min])
    sph = ds.sphere(center,scale)


    thtr = looper.tr
    my_frame = np.where( thtr.frames     == looper.current_frame)[0][0]
    ms = trackage.mini_scrubber(thtr,snapshot.core_id, do_velocity=True)
    if ms.r.shape[0] < 5:
        return -1
    tmap=rainbow_map(ms.ntimes)
    density = thtr.c([snapshot.core_id],'density')
    this_r = ms.r[:,my_frame]
    this_r[ this_r < 1./2048] = 1./2048
    plt.scatter( this_r, density[:,my_frame],c=[tmap(my_frame)]*ms.r.shape[0])

    popt, pcov = curve_fit(powerlaw, this_r, np.log10(density[:,my_frame]), p0=[1,1,-2])
    fit_rho0, fit_r0, fit_alpha = popt
    plt.plot( this_r, 10**powerlaw(this_r, fit_rho0, fit_r0, fit_alpha),c='k')
    looper.alpha_dict_part[looper.current_frame][snapshot.core_id] = popt

    if 1:
        override_bins={}
        override_bins['radius']=np.logspace( np.log10(1./2048), np.log10(0.2),8)
        prof = yt.create_profile(sph,fields[0],fields[1] ,weight_field=weight_field, override_bins=override_bins)#,accumulation=accumulation,
                            #fractional=fractional, n_bins=n_bins, extrema=extrema)
        looper.last_prof=prof
        the_x = 0.5*(prof.x_bins[1:]+prof.x_bins[0:-1])
        the_y = prof[fields[1]]
        ok = the_y > 0
        the_x = the_x[ok]
        the_y = the_y[ok]
        x_units=the_x.units
        y_units=the_y.units
        plt.plot(the_x,the_y,label="n%04d"%looper.current_frame)


        popt, pcov = curve_fit(powerlaw, the_x, np.log10(the_y), p0=[1,1,-2])
        fit_rho0, fit_r0, fit_alpha = popt
        plt.plot( the_x, 10**powerlaw(the_x, fit_rho0, fit_r0, fit_alpha),c='g')
        looper.alpha_dict_prof[looper.current_frame][snapshot.core_id] = popt


        #all_xbins.append(the_x)
        #all_profiles.append(the_y)
        scaledict={True:'log',False:'linear','log':'log','linear':'linear'}
        plt.xscale(scaledict[scales[0]]); plt.yscale(scaledict[scales[1]])
        plt.xlabel(r'%s $%s$'%(fields[0],x_units)); plt.ylabel(r'%s $%s$'%(fields[1],y_units))
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim( ms.r.min(), ms.r.max())
        plt.ylim( density.min(), density.max())
    #plt.legend(loc=0)
    profname = '%s_prof_c%04d_n%04d_%s_%s_%s.png'%(looper.out_prefix, snapshot.core_id, looper.current_frame, fields[0], fields[1], weight_field)
    plt.savefig(profname)
    print(profname)
            
@looper.core_loop
def core_hist(looper,snapshot, field='density', axis_list=[0,1,2], color='r',force_log=None,linthresh=100):
    if snapshot.R_centroid is None:
        snapshot.get_all_properties()
    ds = snapshot.get_ds()
    for ax in axis_list:
        center = ds.arr(snapshot.R_centroid,'code_length')
        Rmax = snapshot.R_mag.max()
        scale_min = ds.arr(0.05,'code_length')
        scale = max([Rmax, scale_min])
        sph = ds.sphere(center,scale)
        prof = yt.create_profile(sph,field,'cell_volume')

from scipy.optimize import curve_fit
def powerlaw(r,rho0, r0, alpha):
    return alpha*np.log10(r/r0) + np.log10(rho0)
@looper.core_loop
def core_radial(looper,snapshot, fields=[],weight_field=None,scales=['log','log'] ):

    if "alpha_dict_part" not in looper.__dict__:
        looper.alpha_dict_part={}
        looper.alpha_dict_prof={}
    if looper.current_frame not in looper.alpha_dict_part:
        looper.alpha_dict_part[looper.current_frame]={}
        looper.alpha_dict_prof[looper.current_frame]={}

    plt.clf()
    if snapshot.R_centroid is None:
        snapshot.get_all_properties()
    ds = snapshot.get_ds()
    center = ds.arr(snapshot.R_centroid,'code_length')
    Rmax = snapshot.R_mag.max()
    scale_min = ds.arr(0.05,'code_length')
    scale = max([Rmax, scale_min])
    sph = ds.sphere(center,scale)


    thtr = looper.tr
    my_frame = np.where( thtr.frames     == looper.current_frame)[0][0]
    ms = trackage.mini_scrubber(thtr,snapshot.core_id, do_velocity=True)
    if ms.r.shape[0] < 5:
        return -1
    tmap=rainbow_map(ms.ntimes)
    density = thtr.c([snapshot.core_id],'density')
    this_r = ms.r[:,my_frame]
    this_r[ this_r < 1./2048] = 1./2048
    plt.scatter( this_r, density[:,my_frame],c=[tmap(my_frame)]*ms.r.shape[0])

    popt, pcov = curve_fit(powerlaw, this_r, np.log10(density[:,my_frame]), p0=[1,1,-2])
    fit_rho0, fit_r0, fit_alpha = popt
    plt.plot( this_r, 10**powerlaw(this_r, fit_rho0, fit_r0, fit_alpha),c='k')
    looper.alpha_dict_part[looper.current_frame][snapshot.core_id] = popt

    if 1:
        override_bins={}
        override_bins['radius']=np.logspace( np.log10(1./2048), np.log10(0.2),8)
        prof = yt.create_profile(sph,fields[0],fields[1] ,weight_field=weight_field, override_bins=override_bins)#,accumulation=accumulation,
                            #fractional=fractional, n_bins=n_bins, extrema=extrema)
        looper.last_prof=prof
        the_x = 0.5*(prof.x_bins[1:]+prof.x_bins[0:-1])
        the_y = prof[fields[1]]
        ok = the_y > 0
        the_x = the_x[ok]
        the_y = the_y[ok]
        x_units=the_x.units
        y_units=the_y.units
        plt.plot(the_x,the_y,label="n%04d"%looper.current_frame)


        popt, pcov = curve_fit(powerlaw, the_x, np.log10(the_y), p0=[1,1,-2])
        fit_rho0, fit_r0, fit_alpha = popt
        plt.plot( the_x, 10**powerlaw(the_x, fit_rho0, fit_r0, fit_alpha),c='g')
        looper.alpha_dict_prof[looper.current_frame][snapshot.core_id] = popt


        #all_xbins.append(the_x)
        #all_profiles.append(the_y)
        scaledict={True:'log',False:'linear','log':'log','linear':'linear'}
        plt.xscale(scaledict[scales[0]]); plt.yscale(scaledict[scales[1]])
        plt.xlabel(r'%s $%s$'%(fields[0],x_units)); plt.ylabel(r'%s $%s$'%(fields[1],y_units))
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim( ms.r.min(), ms.r.max())
        plt.ylim( density.min(), density.max())
    #plt.legend(loc=0)
    profname = '%s_prof_c%04d_n%04d_%s_%s_%s.png'%(looper.out_prefix, snapshot.core_id, looper.current_frame, fields[0], fields[1], weight_field)
    plt.savefig(profname)
    print(profname)

@looper.frame_loop
def slice_raster(self, axis_list=[0,1,2], field='densty',color_dict={}, zlim=None):
    """Full projections, with particles plotted.  Also plots the radius and number for each core.
    Colors can be passed in through a dictionary, indexed by core_id"""
    for axis in axis_list:
        for nslice, zcoord in enumerate( np.arange(0,1,1./128)):
            #proj = self.ds.slice(axis,zcoord, 'c',field)
            #pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
            center = [0.5,0.5,0.5]
            center[axis]=zcoord
            pw = yt.SlicePlot(self.ds, axis, fields=[field],center=center)
            pw.set_cmap(field,'gray')
            pw.set_zlim(field,-1e3,1e3)
            radius_from_core = []
            core_label = ""
            for nc,core_number in enumerate(self.target_indices):
                this_snapshot = self.make_snapshot(self.current_frame,core_number)
                if this_snapshot.R_centroid is None:
                    this_snapshot.get_all_properties()
                center = this_snapshot.R_centroid
                color = color_dict.get(core_number,'r')
                #core_label += "c%04d_"%core_number
                core_label = 'all'


                if np.abs(this_snapshot.R_centroid[axis].v - zcoord) < this_snapshot.R_mag.max().v:
                    pw.annotate_text(cbump(center),
                                     "%d"%core_number,text_args={'color':color}, 
                                     inset_box_args={'visible':False},
                                     coord_system='data')
                    pw.annotate_select_particles(1.0, col=color, indices=self.target_indices[core_number])
            outname = '%s_t0_full_particles_n%04d_sl%04d_%s'%(self.out_prefix,self.current_frame, nslice,core_label)
            print( pw.save(outname))

@looper.core_loop
def core_proj_follow_tweaked(looper,snapshot, field='density', axis_list=[0,1,2], color='r',force_log=None,linthresh=100):
    if snapshot.R_centroid is None:
        snapshot.get_all_properties()
    ds = snapshot.get_ds()
    for ax in axis_list:
        center = ds.arr(snapshot.R_centroid,'code_length')
        Rmax = snapshot.R_mag.max()
        scale_min = ds.arr(0.05,'code_length')
        scale = max([Rmax, scale_min])
        sph = ds.sphere(center,scale)
        proj = ds.proj(field,ax,center=center, data_source = sph) 
        pw = proj.to_pw(center = center,width=(1.0,'code_length'), origin='domain')
        pw.zoom(1./(2*scale.v))
        #pw.zoom(1./0.05)
        pw.set_cmap(field,'gray')
        if force_log is not None:
            pw.set_log(field,force_log,linthresh=linthresh)
        pw.annotate_sphere(center,Rmax, circle_args={'color':color} ) #R_mag.max())
        pw.annotate_text(center,
                         "%d"%snapshot.core_id,text_args={'color':color}, 
                         inset_box_args={'visible':False},
                         coord_system='data')
        thtr = looper.tr
        my_frame = np.where( thtr.frames     == looper.current_frame)[0][0]
        ms = trackage.mini_scrubber(thtr,snapshot.core_id, do_velocity=True)
        pw.annotate_text(center,'X',text_args={'color':'r'})
        other_center = ds.arr([ms.mean_x[my_frame],ms.mean_y[my_frame],ms.mean_z[my_frame]],'code_length')
        pw.annotate_text(center,'O',text_args={'color':'r'})

        pw.annotate_select_particles4(1.0, col=color, indices=snapshot.target_indices)
        pw.annotate_grids()
        if 0:
            pw.set_log('velocity_magnitude',True)
            pw.set_zlim('velocity_magnitude',.1,10)
            pw.set_log('density',True)
            pw.set_zlim('density',.01,10)
        outname = "%s_c%04d_n%04d_centered"%(looper.out_prefix,snapshot.core_id,snapshot.frame)
        print(pw.save(outname))
        fig,ax=plt.subplots(1,1)
        ax0=ax;#[0];ax1=ax[1]
        ax0.scatter( sph['radius'],sph['velocity_magnitude'],c='k')

        my_frame = np.where( looper.tr.frames == looper.current_frame)
        thtr = looper.tr
        ms = trackage.mini_scrubber(thtr,snapshot.core_id, do_velocity=True)
        ax0.scatter( ms.r[:,my_frame], ms.raw_v2[:,my_frame]**0.5,c='r')
        fig.savefig('plots_to_sort/scatter_c%04d_n%04d.png'%(snapshot.core_id,looper.current_frame))
        plt.close(fig)
    return pw
    
