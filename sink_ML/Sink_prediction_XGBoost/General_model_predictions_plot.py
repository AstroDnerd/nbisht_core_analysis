# standard system modules
import os, sys
os.environ["PATH"] += os.pathsep + "/home/nbisht/myapps/bin/"
import h5py 
import pandas as pd
import pickle
# standard module for array manipulation
import numpy as np

# standard module for high-quality plots
from PIL import Image
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
mp.rcParams.update(mp.rcParamsDefault)


# set a seed to ensure reproducibility
seed = 128
rnd  = np.random.RandomState(seed)

def mae_modded(y_true, y_pred):
    mae = np.array([0.,0.,0.])
    mae_diff = np.abs(y_true - y_pred)
    mae_add = 1 - mae_diff
    stacked = np.stack([mae_diff, mae_add], axis=2)
    mae = stacked.min(axis=2).mean(axis=0)
    return mae

def r2_modded(y_true, y_pred):
    true_mean = np.mean(y_true,axis=0)
    mae_diff = np.abs(y_true - y_pred)
    mae_add = 1 - mae_diff
    stacked = np.stack([mae_diff, mae_add], axis=2)
    r2_residual = stacked.min(axis=2)**2
    mae_diff = np.abs(y_true - true_mean)
    mae_add = 1 - mae_diff
    stacked = np.stack([mae_diff, mae_add], axis=2)
    r2_total = stacked.min(axis=2)**2

    return 1-r2_residual.sum(axis=0)/r2_total.sum(axis=0)

with open('/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_Core_predictions.pickle', 'rb') as handle:
    results_dic_core = pickle.load(handle)

with open('/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_NonCore_predictions.pickle', 'rb') as handle:
    results_dic_noncore = pickle.load(handle)

with open('/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_Combined_predictions.pickle', 'rb') as handle:
    results_dic_combined = pickle.load(handle)


#{'Model_1':{'ytrue':[], 'ypred':[]}, 'Model_2':{'ytrue':[], 'ypred':[]}, 'Model_3':{'ytrue':[], 'ypred':[]}}

TARGET = ['X', 'Y', 'Z']
prediction_type = [results_dic_core, results_dic_noncore, results_dic_combined]
model_names = ['Model_1', 'Model_2', 'Model_3']
fig = plt.figure(figsize=(27, 9))
for i in range(3): #loop over models
    for j in range(3): #loop over prediction types
        y_pred = prediction_type[j][str(model_names[i])+'_ypred'].to_numpy()
        y_true = prediction_type[j]['ytrue'].to_numpy()
        mae = mae_modded(y_true, y_pred)
        r2 = r2_modded(y_true, y_pred)
        for k in range(3): #loop over prediction axes
            ax = plt.subplot(3, 9, i*9+j*3+k+1)
            ax.scatter(y_true[:,k], y_pred[:,k], c='crimson', s=1e-4)
            ax.plot([0, 0], [1, 1], 'b-')
            ax.set_title(f'{TARGET[k]}: MAE = {mae[k]:.6f}, R2 = {r2[k]:.4f}', fontsize=12)
            ax.set_xlim([-0.1,1.1])
            ax.set_ylim([-0.1,1.1])
fig.tight_layout()
plt.savefig('./sink_ML/Sink_prediction_XGBoost/Best_Models_General_prediction.png')
plt.show()