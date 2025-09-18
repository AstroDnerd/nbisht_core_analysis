

all_cores = leaf_indices.keys()
n_particles = np.array([len(leaf_indices[ic]) for ic in all_cores])
plt.clf()
plt.plot(all_cores,n_particles,marker='*')
plt.ylabel('N_particles')
plt.xlabel('core id')
plt.savefig('n_cores.pdf')


