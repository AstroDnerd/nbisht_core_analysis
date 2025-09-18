import imageio
import sys
import os

def animator(path_to_plots,name_of_gif):
	images = []
	filenames = list(sorted(os.listdir(path_to_plots)))

	for filename in filenames:
	    images.append(imageio.imread(path_to_plots+'/'+filename))
	imageio.mimsave(path_to_plots+'/../'+name_of_gif+'.gif', images)