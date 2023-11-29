import bigfish.stack as stack
import bigfish.plot as plot
import bigfish.detection as detection
import numpy as np
import matplotlib.pyplot as plt

path = '/media/floricslimani/SSD 4To/SSD_floricslimani/4_centrosome/input/r06c05-HeLa-c_inh216h_KIF1Cf06-EGFP-sk1fk1fl1.tiff'
show= False

voxel_size = (300, 103, 103)
spot_size = (800,800,800)


egfp = stack.read_image(path)
log = stack.log_filter(egfp, sigma= 2)
max_proj = stack.maximum_projection(egfp)
mean_proj = stack.mean_projection(egfp)

plot.plot_images([max_proj, stack.mean_projection(log), stack.maximum_projection(log)],titles= ['max', 'log_mean', 'log_max'], rescale= True)