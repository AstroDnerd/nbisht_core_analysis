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

#plt.style.use(['science','scatter','grid'])

# set a seed to ensure reproducibility
seed = 128
rnd  = np.random.RandomState(seed)
FRAME_DIFF = 30


TARGET = ['X', 'Y', 'Z']
prediction_type_name = ['Core', 'NonCore', 'Combined']
model_names = ['Model_1', 'Model_2', 'Model_3']


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

from matplotlib.lines import Line2D
from matplotlib.patches import Patch

def paper_plot_mae_r2_phase_plot(X_test, ypred, ytrue, append=''):
    axes_colors = ['tomato', 'limegreen', 'dodgerblue']
    axes_colors_pred = ['grey']
    fig, ax_arr = plt.subplots(1, 3, figsize=(15, 4))
    y_true = results_dic_core['ytrue'].to_numpy()
    y_pred = results_dic_core['Model_1_ypred'].to_numpy()
    x_test = results_dic_core['xtest'][['X_i', 'Y_i', 'Z_i']].to_numpy()
    mae = mae_modded(y_true, y_pred)
    r2 = r2_modded(y_true, y_pred)
    for k in range(3):
        ax = ax_arr[k]
        ax.axline([0, 0], [1, 1], color='black', linestyle='--', lw=0.5)
        ax.scatter(y_true[:,k], x_test[:,k], c=axes_colors[k], s=1e-4, alpha = 0.5, label='Truth', zorder = 0)
        ax.scatter(y_true[:,k], y_pred[:,k], c=axes_colors_pred, s=1e-4, alpha = 1, label='Prediction', zorder=10)
        ax.text(0.5, 1.01, f'{TARGET[k]}: \n MAE = {mae[k]:.6f} \n R2 = {r2[k]:.4f}', fontsize=8, horizontalalignment='center', verticalalignment='top')
        ax.set_aspect('equal')
        ax.set_ylabel('Initial', fontsize=12)
        #ax.set_ylabel(f'Core {model_names[0]} Prediction', fontsize=12)
        ax.set_xlabel('Final', fontsize=12)

    legend_elements = [Line2D([0], [0], marker='o', color=axes_colors[0], label=r'$X$', markersize=1),
                        Line2D([0], [0], marker='o', color=axes_colors[1], label=r'$Y$', markersize=1),
                        Line2D([0], [0], marker='o', color=axes_colors[2], label=r'$Z$', markersize=1),
                        Line2D([0], [0], marker='o', color='0.5', label='Prediction', markersize=1)
                        ]

    fig.legend(handles=legend_elements, loc='upper center', fontsize=10, ncol=4)
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    plt.savefig('./sink_ML/Sink_prediction_XGBoost/Best_Models_General_prediction'+append+'_initial.png', dpi=300, bbox_inches='tight')
    plt.close()

if 0:
    with open('/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/u502_framewise_predictions_using_2_sims.pickle', 'rb') as handle:
        results_dic_core = pickle.load(handle)
    col = 0
    row = 0
    prediction_type = [results_dic_core, results_dic_core, results_dic_core]
    ypred = prediction_type[col][str(model_names[row])+'_ypred'].reset_index(drop=True)
    ytrue = prediction_type[col]['ytrue'].reset_index(drop=True)
    X_test = prediction_type[col]['xtest'].reset_index(drop=True)
    paper_plot_mae_r2_phase_plot(X_test, ypred, ytrue, append='_u502_2sims')

def plot_prediction_framewise(X_test, ypred, ytrue, append=''):
    unique_frames = np.unique(X_test['Initial_Frame'])
    if len(unique_frames)%2 == 0:
        num_cols = len(unique_frames)//2
    else:
        num_cols = len(unique_frames)//2+1
    fig = plt.figure(figsize=(4.5*num_cols, 12))
    for frame_num_index in range(0,len(unique_frames),2):
        frame_val = unique_frames[frame_num_index]
        print(frame_num_index,frame_val)
        X_test_frame = X_test[X_test['Initial_Frame'] == frame_val]
        ypred_frame = ypred.loc[X_test_frame.index]
        ytrue_frame = ytrue.loc[X_test_frame.index]
        ax = plt.subplot2grid((3,num_cols), (0,frame_num_index//2))
        ax.set_title(f'Initial Frame: {frame_val}', fontsize=18)
        ax.set_xlim(0, 1)
        ax.set_xlabel('X', fontsize=18)
        ax.set_ylim(0, 1)
        ax.set_ylabel('Y', fontsize=18)
        ax.scatter(X_test_frame['X_i'], X_test_frame['Y_i'], s=1e-3, color='coral', alpha = 0.5)

        ax = plt.subplot2grid((3,num_cols), (1,frame_num_index//2))
        ax.set_title(f'True Frame: {frame_val+FRAME_DIFF}', fontsize=18)
        ax.set_xlim(0, 1)
        ax.set_xlabel('X', fontsize=18)
        ax.set_ylim(0, 1)
        ax.set_ylabel('Y', fontsize=18)
        ax.scatter(ytrue_frame['X_f'], ytrue_frame['Y_f'], s=1e-3, color='cornflowerblue', alpha = 0.5)

        ax = plt.subplot2grid((3,num_cols), (2,frame_num_index//2))
        ax.set_title(f'Predicted Frame: {frame_val+FRAME_DIFF}', fontsize=18)
        ax.set_xlim(0, 1)
        ax.set_xlabel('X', fontsize=18)
        ax.set_ylim(0, 1)
        ax.set_ylabel('Y', fontsize=18)
        ax.scatter(ypred_frame['X_f'], ypred_frame['Y_f'], s=1e-3, color='royalblue', alpha = 0.5)
        
    fig.tight_layout()
    plt.savefig('./sink_ML/Sink_prediction_XGBoost/Best_Model1_Framewise_prediction_with_initialframe_'+append+'.png', dpi=300, bbox_inches='tight')
    plt.show()


def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def set_axis_style(ax, labels):
    ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Sample name')

def plot_framewise_l2norm(X_test, ypred, ytrue, append=''):
    unique_frames = np.unique(X_test['Initial_Frame'])
    predicted_dataset = []
    predicted_as_orig_dataset = []
    predicted_frame = []
    for frame_num_index in range(0,len(unique_frames)):
        frame_val = unique_frames[frame_num_index]
        predicted_frame.append(frame_val+FRAME_DIFF)
        X_test_frame = X_test[X_test['Initial_Frame'] == frame_val]
        X_initial = X_test_frame[['X_i', 'Y_i', 'Z_i']]
        X_initial = X_initial.loc[X_test_frame.index]
        X_initial.rename(columns={"X_i": "X_f", "Y_i": "Y_f", 'Z_i':'Z_f'}, inplace=True)
        ypred_frame = ypred.loc[X_test_frame.index]
        ytrue_frame = ytrue.loc[X_test_frame.index]

        diff_orig = np.abs(ytrue_frame-X_initial)
        diff = np.abs(ytrue_frame-ypred_frame)
        for df in [diff_orig, diff]:
            df['X_f'] = np.where(df['X_f']>=0.5, 1-df['X_f'], df['X_f'])
            df['Y_f'] = np.where(df['Y_f']>=0.5, 1-df['Y_f'], df['Y_f'])
            df['Z_f'] = np.where(df['Z_f']>=0.5, 1-df['Z_f'], df['Z_f'])
        diff = diff.to_numpy()
        diff_orig = diff_orig.to_numpy()
        predicted_dataset.append(np.sort(np.linalg.norm(diff, axis=1)))
        predicted_as_orig_dataset.append(np.sort(np.linalg.norm(diff_orig, axis=1)))
    data = np.array(predicted_dataset).T
    data_orig = np.array(predicted_as_orig_dataset).T
    print(data.shape)
    print(len(predicted_frame))
    fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(24, 9))
    parts = ax1.violinplot(
            data, showmeans=False, showmedians=False,
            showextrema=False)

    for pc in parts['bodies']:
        pc.set_facecolor('plum')
        pc.set_edgecolor('indigo')
        pc.set_alpha(1)

    quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=0)
    medians_orig= np.percentile(data_orig,  50, axis=0)
    whiskers = np.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(data, quartile1, quartile3)])
    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]

    inds = np.arange(1, len(medians) + 1)
    #prediction norm
    ax1.scatter(inds, medians, marker='o', color='white', s=30, zorder=3)
    ax1.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
    ax1.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

    #original as predicted
    ax1.scatter(inds, medians_orig, marker='x', color='red', s=50, zorder=3)
    print(medians_orig)

    ax1.axhline(0, color='black', linestyle='--', lw=2)
    ax1.set_ylabel('Euclidean Distance', fontsize=18)
    # set style for the axes
    for ax in [ax1]:
        set_axis_style(ax, predicted_frame)
    
    ax1.set_xlabel('Frame', fontsize=18)


    plt.savefig('./sink_ML/Sink_prediction_XGBoost/Best_Model1_Framewise_Euclidean_Distance'+append+'.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


#{'Model_1':{'ytrue':[], 'ypred':[]}, 'Model_2':{'ytrue':[], 'ypred':[]}, 'Model_3':{'ytrue':[], 'ypred':[]}}

if 0:
    with open('/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_Core_framewise_predictions.pickle', 'rb') as handle:
        results_dic_core = pickle.load(handle)
    with open('/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_NonCore_predictions.pickle', 'rb') as handle:
        results_dic_noncore = pickle.load(handle)
    with open('/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_Combined_predictions.pickle', 'rb') as handle:
        results_dic_combined = pickle.load(handle)
    col = 0
    row = 0
    prediction_type = [results_dic_core, results_dic_core, results_dic_core]
    ypred = prediction_type[col][str(model_names[row])+'_ypred'].reset_index(drop=True)
    ytrue = prediction_type[col]['ytrue'].reset_index(drop=True)
    X_test = prediction_type[col]['xtest'].reset_index(drop=True)
    plot_framewise_l2norm(X_test, ypred, ytrue)

if 1:
    with open('/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/u502_framewise_predictions_using_2_sims.pickle', 'rb') as handle:
        results_dic_core = pickle.load(handle)
    col = 0
    row = 0
    prediction_type = [results_dic_core, results_dic_core, results_dic_core]
    ypred = prediction_type[col][str(model_names[row])+'_ypred'][['X_f', 'Y_f', 'Z_f']].reset_index(drop=True)
    ytrue = prediction_type[col]['ytrue'][['X_f', 'Y_f', 'Z_f']].reset_index(drop=True)
    X_test = prediction_type[col]['xtest'].reset_index(drop=True)
    plot_prediction_framewise(X_test, ypred, ytrue,append='_u502_2sims')
    plot_framewise_l2norm(X_test, ypred, ytrue, append='_u502_2sims')



