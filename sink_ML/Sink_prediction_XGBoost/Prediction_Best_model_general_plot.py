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
import matplotlib.gridspec as gridspec
import scienceplots

#plt.style.use(['science','scatter','grid'])

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

def plot_all_true_pred_panels():
    with open('/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_Core_predictions.pickle', 'rb') as handle:
        results_dic_core = pickle.load(handle)

    with open('/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_NonCore_predictions.pickle', 'rb') as handle:
        results_dic_noncore = pickle.load(handle)

    with open('/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_Combined_predictions.pickle', 'rb') as handle:
        results_dic_combined = pickle.load(handle)


    #{'Model_1':{'ytrue':[], 'ypred':[]}, 'Model_2':{'ytrue':[], 'ypred':[]}, 'Model_3':{'ytrue':[], 'ypred':[]}}

    TARGET = ['X', 'Y', 'Z']
    prediction_type_name = ['Core', 'NonCore', 'Combined']
    prediction_type = [results_dic_core, results_dic_noncore, results_dic_combined]
    model_names = ['Model_1', 'Model_3']
    axes_colors = ['tomato', 'limegreen', 'dodgerblue']
    fig = plt.figure(figsize=(24, 6))
    outer = gridspec.GridSpec(nrows=2, ncols=3, hspace=0.2)
    mae_dic = np.zeros((2,3))
    r2_dic = np.zeros((2,3))
    for row in range(2):
        for col in range(3):
            inner = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=outer[row, col], wspace=0)
            y_pred = prediction_type[col][str(model_names[row])+'_ypred'].to_numpy()
            y_true = prediction_type[col]['ytrue'].to_numpy()
            mae = mae_modded(y_true, y_pred)
            r2 = r2_modded(y_true, y_pred)
            mae_dic[row, col] = mae.mean()
            r2_dic[row, col] = r2.mean()
            axs = [plt.subplot(cell) for cell in inner]
            for k in range(len(axs)):
                ax = axs[k]
                if k==1 and row==0:
                    ax.set_title(f'{prediction_type_name[col]}', fontsize=18)
                if k==1 and row==1:
                    ax.set_xlabel('Truth', fontsize=12)
                if k>0:
                    ax.set_yticklabels([])
                    ax.set_yticks([])
                if row==0:
                    ax.set_xticklabels([])
                ax.axline([0, 0], [1, 1], color='black', linestyle='--', lw=0.5)
                ax.scatter(y_true[:,k], y_pred[:,k], c=axes_colors[k], s=1e-4, alpha = 0.5)
                ax.text(0.5, 1.08, f'{TARGET[k]}: \n MAE = {mae[k]:.6f} \n R2 = {r2[k]:.4f}', fontsize=10, horizontalalignment='center', verticalalignment='top')
                ax.set_xlim([-0.1,1.1])
                ax.set_ylim([-0.1,1.1])
                ax.set_aspect('equal')
                if col==0 and k==0:
                    ax.set_ylabel(f'{model_names[row]} Prediction', fontsize=12)

    plt.show()
    plt.savefig('./sink_ML/Sink_prediction_XGBoost/Best_Models_General_prediction.png', dpi=300, bbox_inches='tight')
    plt.close()
    return mae_dic, r2_dic

#mae_dic, r2_dic = plot_all_true_pred_panels()
#print(mae_dic)
#print(r2_dic)
mae_dic = [[0.00680006, 0.00799723, 0.00638501], #Model 1
            [0.00674917, 0.00797059, 0.00632831]] #Model 3

r2_dic = [[0.99882138, 0.99853539, 0.99895557],
            [0.9988365,  0.99853752, 0.99897092]]

def plot_mae_r2_phase_plot(mae_dic,r2_dic):
    with open('/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_Core_predictions.pickle', 'rb') as handle:
        results_dic_core = pickle.load(handle)

    TARGET = ['X', 'Y', 'Z']
    prediction_type_name = ['Core', 'NonCore', 'Combined']
    model_names = ['Model 1', 'Model 3']
    axes_colors = ['tomato', 'limegreen', 'dodgerblue']
    fig, ax_arr = plt.subplots(3, 1, figsize=(6, 18), sharey=True)
    y_true = results_dic_core['ytrue'].to_numpy()
    y_pred = results_dic_core['Model_1_ypred'].to_numpy()
    mae = mae_modded(y_true, y_pred)
    r2 = r2_modded(y_true, y_pred)
    ax_arr[0].set_title(f'{prediction_type_name[0]}', fontsize=18)
    for k in range(3):
        ax = ax_arr[k]
        ax.axline([0, 0], [1, 1], color='black', linestyle='--', lw=0.5)
        ax.scatter(y_true[:,k], y_pred[:,k], c=axes_colors[k], s=1e-4, alpha = 0.5)
        ax.text(0.5, 1.08, f'{TARGET[k]}: \n MAE = {mae[k]:.6f} \n R2 = {r2[k]:.4f}', fontsize=10, horizontalalignment='center', verticalalignment='top')
        ax.set_xlim([-0.1,1.1])
        ax.set_ylim([-0.1,1.1])
        ax.set_aspect('equal')
    
    ax_arr[1].set_ylabel(f'{model_names[0]} Prediction', fontsize=12)
    ax_arr[2].set_xlabel('Truth', fontsize=12)

    fig.subplots_adjust(wspace=0, hspace=0)
    #fig.tight_layout()
    plt.show()
    plt.savefig('./sink_ML/Sink_prediction_XGBoost/Best_Model1_General_prediction.png', dpi=300, bbox_inches='tight')
    plt.close()

    model_colors = ['slateblue', 'orchid']
    scatter_type = ['o', '+', '*']
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    for i in range(2):
        for j in range(3):
            ax.scatter(mae_dic[i][j], r2_dic[i][j], label=prediction_type_name[j] + ' ' +model_names[i], color=model_colors[i], marker=scatter_type[j])
    
    ax.set_xlabel(r'$MAE_{modified}$', fontsize=12)
    ax.set_ylabel(r'$R2_{modified}$', fontsize=12)
    ax.axhline(1, color='black', linestyle='--', lw=0.5)
    ax.legend()
    plt.show()
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    plt.savefig('./sink_ML/Sink_prediction_XGBoost/Best_Models_General_prediction_Phase_Plot.png', dpi=300, bbox_inches='tight')
    plt.close()

def paper_plot_mae_r2_phase_plot(mae_dic,r2_dic):
    with open('/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_Core_predictions.pickle', 'rb') as handle:
        results_dic_core = pickle.load(handle)

    TARGET = ['X', 'Y', 'Z']
    prediction_type_name = ['Core', 'NonCore', 'Combined']
    model_names = ['Model 1', 'Model 3']
    axes_colors = ['tomato', 'limegreen', 'dodgerblue']
    fig, ax_arr = plt.subplots(1, 2, figsize=(15, 6))
    y_true = results_dic_core['ytrue'].to_numpy()
    y_pred = results_dic_core['Model_1_ypred'].to_numpy()
    mae = mae_modded(y_true, y_pred)
    r2 = r2_modded(y_true, y_pred)
    k=0
    ax = ax_arr[k]
    ax.axline([0, 0], [1, 0], color='black', linestyle='-.', lw=0.5)
    diff_arr = np.abs(y_true[:,k]- y_pred[:,k])
    add_arr = 1 - diff_arr
    stacked = np.stack([diff_arr, add_arr])
    min_val = stacked.min(axis=0)
    ax.scatter(y_true[:,k], y_true[:,k]- y_pred[:,k], c=axes_colors[k], s=1e-4, alpha = 0.5)
    bins = np.linspace(0, 1, 50)
    digitized = np.digitize(y_true[:,k], bins)
    bin_means = [min_val[digitized == i].mean() for i in range(1, len(bins))]
    ax.plot(bins[1:], bin_means, c = 'slateblue', lw = 2)
    ax.text(0.5, 1.08, f'{TARGET[k]}: \n MAE = {mae[k]:.6f} \n R2 = {r2[k]:.4f}', fontsize=10, horizontalalignment='center', verticalalignment='top')
    ax.set_xlim([-0.05,1.05])
    ax.set_ylim([-1.1,1.1])
    #ax.set_aspect('equal')
    ax.set_ylabel(f'Core {model_names[0]} Prediction - Truth', fontsize=12)
    ax.set_xlabel('Truth', fontsize=12)

    ax = ax_arr[1]
    model_colors = ['slateblue', 'orchid']
    scatter_type = ['o', '+', '*']
    for i in range(2):
        for j in range(3):
            ax.scatter(mae_dic[i][j], r2_dic[i][j], label=prediction_type_name[j] + ' ' +model_names[i], color=model_colors[i], marker=scatter_type[j])
    
    ax.set_xlabel(r'$MAE_{modified}$', fontsize=12)
    ax.set_ylabel(r'$R2_{modified}$', fontsize=12)
    ax.axhline(1, color='black', linestyle='--', lw=0.5)
    ax.legend()
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    plt.savefig('./sink_ML/Sink_prediction_XGBoost/Best_Models_General_prediction_error.png', dpi=300, bbox_inches='tight')
    plt.close()
paper_plot_mae_r2_phase_plot(mae_dic,r2_dic)