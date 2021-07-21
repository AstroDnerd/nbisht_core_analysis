from starter2 import *


def shift_snaps(loop):
    for frame in loop.snaps:
        for core_id in loop.snaps[frame]:
            snap = loop.snaps[frame][core_id]
            snap.ds = loop.load(frame)
            tr = loop.tr
            ms = trackage.mini_scrubber( tr, core_id, do_velocity=False)
            ms_index = np.where( tr.frames == frame )[0]

            particle_ids = loop.tr.c([core_id],'particle_id')
            if hasattr(particle_ids,'v'):
                particle_ids = particle_ids.v
            total_particle_index_error = np.abs( particle_ids  - snap.ind.v).sum()
            if total_particle_index_error > 0:
                t1 = np.abs( particle_ids - np.sort(particle_ids)).sum()
                t2 = np.abs( snap.ind.v - np.sort(snap.ind.v)).sum()
                print("Particle ids are sorted if zero:", t1)
                print("Snap.ind are sorted if zero:", t2)
                raise("SORT ERROR")
            pos2 = np.concatenate([ms.this_x[:,ms_index], ms.this_y[:,ms_index], ms.this_z[:,ms_index]], axis=1)

            maxshift = np.abs( snap.pos.v - pos2).max() 
            print( "MAX POS SHIFT", frame, core_id, maxshift)
            if maxshift > 0.1:
                shift =  pos2 - snap.pos.v
                shift_amount = snap.ds.arr(np.sign(shift), 'code_length')
                to_shift = np.abs( shift) > 0.25
                #if to_shift.sum():
                #    pdb.set_trace()
                print("TO_SHIFT", to_shift)
                snap.pos[ to_shift ] += shift_amount[ to_shift ]
                print("SHIFTED TO",snap.pos[to_shift], shift_amount[to_shift])
            density = snap.field_values['density']
            volume = snap.field_values['cell_volume']
            m = (density*volume).sum()
            centroid_tmp =  np.array([(snap.pos[:,dim]*density*volume).sum()/m for dim in range(3)])
            snap.R_centroid = snap.ds.arr(centroid_tmp,'cm')
            print("SNAP ", snap.R_centroid)
            print(pos2.min())

            #
            # OK for next time, snap.pos and pos2 should all be on the same side of the box
            #
