import pyximport; pyximport.install()
import particle_ops
import particle_grid_mask

def get_current_mask_test(self):
    """get the particle mask that relates particles for this core_id and this frame
    to the particles in the target_indices from the target_frame"""
    these_pids=self.target_indices.astype('int64')
    if type(these_pids) == yt.units.yt_array:
        these_pids = these_pids.v
    data_region = self.get_region(self.frame)
    mask_to_get=np.zeros(these_pids.shape,dtype='int32')
    my_indices = data_region['particle_index'].astype('int64')
    bad_particles=[]
    for ppp in these_pids:
        if ppp not in my_indices:
            print("SHIT particle %d not found"%ppp)
            bad_particles.append(ppp)
    #found_any, mask = particle_ops.mask_particles(these_pids,my_indices,mask_to_get)
    #self.mask = mask
    #return found_any, mask
    return bad_particles
#bad_particles = get_current_mask_test(sn)
grids=[]
for grid in sn.ds.index.grids:
    for ppp in bad_particles:
        if ppp in grid['particle_index']:
            ind = np.where( grid['particle_index']==ppp)[0]
            print("horray found it %d"%ppp)
            print(grid)
            pos=grid['particle_position'][ind]
            print( "Left  %s"%(str(grid.LeftEdge)))
            print( "Right %s"%(str(grid.RightEdge)))
            print( "pos   %s"%(str(pos)))

            print( (pos-grid.LeftEdge)/grid.dds)
            print( (pos-grid.RightEdge)/grid.dds)
            grids.append(grid)
