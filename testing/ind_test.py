
left = np.array([1,1,1])
right = np.array([0,0,0])
for core_id in core_list:
    snapshot = tl.looper1.snaps[frame][core_id]
    this_left = snapshot.pos.min(axis=0)
    this_right = snapshot.pos.max(axis=0)
    left = np.row_stack([this_left,left]).min(axis=0)
    right = np.row_stack([this_right,right]).max(axis=0)


