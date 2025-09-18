# standard system modules
import os, sys
os.environ["PATH"] += os.pathsep + "/home/nbisht/myapps/bin/"
import h5py 
import pandas as pd
import pickle
# standard module for array manipulation
import numpy as np
import xgboost as xgb
# standard module for high-quality plots
from PIL import Image
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
mp.rcParams.update(mp.rcParamsDefault)
import matplotlib.gridspec as gridspec
import scienceplots
from io import BytesIO
import PIL
#plt.style.use(['science','scatter','grid'])

# set a seed to ensure reproducibility
seed = 128
rnd  = np.random.RandomState(seed)


def paper_plot_ensemble_tree(model):
    #g = xgb.to_graphviz(model, num_trees=0, ax=ax,graph_attrs= {"with_stats": True, 'rankdir':'LR'}).render('./plots_to_sort/Ensemble_tree_latex', format='xdot')
    for i in range(150):
        fig,ax = plt.subplots(figsize=(15, 15))
        g = xgb.to_graphviz(model, num_trees=i, ax=ax,graph_attrs= {"with_stats": True})
        s = BytesIO()
        s.write(g.pipe(format="eps"))
        s.seek(0)
        img = np.array(PIL.Image.open(s))
        if img.shape[1] < 3000:
            print(i, img.shape)
            plt.imshow(img)
            plt.axis("off")
            plt.savefig('./plots_to_sort/Ensemble_trees/Ensemble_tree_'+str(i)+'.png', dpi=1200, bbox_inches='tight')
        plt.close()
    


model = xgb.XGBRegressor()
mtype = 'Core'
model.load_model('./sink_ML/Sink_prediction_XGBoost/'+mtype+'_Extended_Framewise_BestModel0_nnmodel.json')
paper_plot_ensemble_tree(model)