
core_cmap = colors.make_core_cmap( new_looper.core_list)
loop_apps.core_proj_multiple(new_looper,core_list=[10],
                             #frame_list=[0],
                             axis_list=[1], color_dict=core_cmap, particles=True,
                             only_sphere=True,zoom=True,
                             center_on_sphere=True, annotate=True,tracker_positions=True, shifted_tracker=False)
