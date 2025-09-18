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
mp.rcParams['agg.path.chunksize'] = 100000
import matplotlib.gridspec as gridspec
import scienceplots

plt.style.use(['science','grid'])

# set a seed to ensure reproducibility
seed = 128
rnd  = np.random.RandomState(seed)
FRAME_DIFF = 30
#Average Time difference between snapshots
DELTAT = 0.0003744333333
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


with open('/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_Core_corewise_predictions.pickle', 'rb') as handle:
    results_dic_core = pickle.load(handle)


#{'Model_1':{'ytrue':[], 'ypred':[]}, 'Model_2':{'ytrue':[], 'ypred':[]}, 'Model_3':{'ytrue':[], 'ypred':[]}}

Test_cores = [229,  15,  59,  58, 203,  34,  96,  10, 106,  17, 238,  98,  23, 134, 150, 148, 188, 105, 135, 191, 259, 252, 225, 220, 284, 262, 230, 154, 157, 156,  91,  71, 171, 269, 129, 162,  
            40,  47,  48, 39,  38, 183, 165,  65, 112, 235, 245, 115, 122, 164,  76, 120, 264, 272, 127, 186,   7]
TARGET = ['X', 'Y', 'Z']
prediction_type_name = ['Core', 'NonCore', 'Combined']
prediction_type = [results_dic_core, results_dic_core, results_dic_core]
model_names = ['Model_1', 'Model_2', 'Model_3']

def plot_corewise_tracks(df_timeseries_test, ypred, ytrue):
    df_timeseries_test_with_pred = df_timeseries_test.copy()
    df_timeseries_test_with_pred['X_f_pred'] = ypred['X_f']
    df_timeseries_test_with_pred['Y_f_pred'] = ypred['Y_f']
    df_timeseries_test_with_pred['Z_f_pred'] = ypred['Z_f']
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    for chosen_core in [229, 96, 238, 245, 47]:
        df_this_core = df_timeseries_test_with_pred[df_timeseries_test_with_pred['Core_id']==chosen_core]
        unique_pids = np.unique(df_this_core['Particle_id'])
        true_tracks = []
        pred_tracks = []
        frames = np.array(df_this_core[df_this_core['Particle_id']==unique_pids[0]]['Initial_Frame'].values)
        frames+=FRAME_DIFF
        all_true_tracks = []
        all_pred_tracks = []
        if len(unique_pids)>1000:
            unique_pids = unique_pids[::2]
        for pid in unique_pids:
            df_this_pid = df_this_core[df_this_core['Particle_id']==pid]
            true_tracks.append(df_this_pid[['X_f', 'Y_f', 'Z_f']].values.T)
            pred_tracks.append(df_this_pid[['X_f_pred', 'Y_f_pred', 'Z_f_pred']].values.T)
            ax.plot(frames*DELTAT, true_tracks[-1][0], 'lightcoral', lw=0.5, alpha = 0.5)
            ax.plot(frames*DELTAT, pred_tracks[-1][0], 'plum', lw=0.5, alpha = 0.5)
            all_true_tracks.append(true_tracks[-1][0])
            all_pred_tracks.append(pred_tracks[-1][0])
        
        median_true = np.median(all_true_tracks, axis=0)
        median_pred = np.median(all_pred_tracks, axis=0)
        ax.plot(frames*DELTAT, median_true, 'firebrick', lw=1, alpha = 1, zorder = 10)
        ax.plot(frames*DELTAT, median_pred, 'rebeccapurple', lw=1, alpha = 1, zorder = 10)
        ax.text(frames[-7]*DELTAT, median_true[-7]+0.01, 'Core '+str(chosen_core), fontsize=5)

    ax.set_xlabel('Time (s)', fontsize=12)
    ax.set_ylabel('X (cm)', fontsize=12)
    ax.set_ylim(0, 1)
    ax.set_xlim(50*DELTAT, 119*DELTAT)
    plt.savefig('./sink_ML/Sink_prediction_XGBoost/Best_Model1_Corewise_tracks.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


col = 0
row = 0
ypred = prediction_type[col][str(model_names[row])+'_ypred'].reset_index(drop=True)
ytrue = prediction_type[col]['ytrue'].reset_index(drop=True)
df_timeseries_test = prediction_type[col]['test_df'].reset_index(drop=True)
plot_corewise_tracks(df_timeseries_test, ypred, ytrue)




