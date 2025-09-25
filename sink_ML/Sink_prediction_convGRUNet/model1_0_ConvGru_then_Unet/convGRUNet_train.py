# standard system modules
import os, sys
import h5py 
import argparse
# standard module for tabular data
import pandas as pd

# standard module for array manipulation
import numpy as np
from itertools import permutations

# standard statistical module
import scipy.stats as st
from scipy import linalg
from scipy.stats import ks_2samp
import skimage as ski

# standard module for high-quality plots
from PIL import Image
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
mp.rcParams.update(mp.rcParamsDefault)

# to plot pixelized images
import imageio.v3 as im

# standard research-level machine learning toolkit from Meta (FKA: FaceBook)
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torch.utils.data as data
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import normalize
import tables
import torchvision

from torchvision.utils import save_image
from model1_0_modules import *

from tqdm import tqdm

# set a seed to ensure reproducibility
seed = 128
rnd  = np.random.RandomState(seed)

DATAFILE  = '/data/nbisht/projects/128/B2/datasets/nb101_ML_dataset_AllData_AutoEnc.h5'
DATASETNAME = '/data/nbisht/projects/128/B2/datasets/Unet/nb101_Images_seq2_3D.hdf5'
MODELFILE = 'nnmodel.dict'

IMAGESIZE = 64
FRAME_DIFF = 30
TEST_PERCENTAGE = 0.01

DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
#DEVICE = torch.device('cpu')

print(f'Available device: {str(DEVICE):4s}', flush=True)

def model_performance(pred_output_arr, act_output_arr):
    # Calculate metrics
    pred_MSE = ski.metrics.mean_squared_error(act_output_arr, pred_output_arr)

    inversed_pred_output_arr = img_inverse_transform(pred_output_arr)
    inversed_act_output_arr = img_inverse_transform(act_output_arr)

    pred_total_density = np.sum(inversed_pred_output_arr)
    act_total_density = np.sum(inversed_act_output_arr)
    density_difference = np.abs(pred_total_density - act_total_density)

    return pred_MSE, density_difference

def plot_comparison(input_arr, pred_output_arr, act_output_arr, label, save_filename='xy_data.png', dont_show = False):
    input_str = 'INPUT IMAGE:\n'
    pred_str = 'PREDICTED IMAGE:\n'
    act_str = 'ACTUAL IMAGE:\n'
    # Calculate metrics
    input_MSE = ski.metrics.mean_squared_error(act_output_arr, input_arr)
    input_str+= 'MSE: '+'{:0.5f}'.format(input_MSE)+'\n'
    pred_MSE = ski.metrics.mean_squared_error(act_output_arr, pred_output_arr)
    pred_str+= 'MSE: '+'{:0.5f}'.format(pred_MSE)+'\n'
    actual_MSE = ski.metrics.mean_squared_error(act_output_arr, act_output_arr)
    act_str+= 'MSE: '+'{:0.5f}'.format(actual_MSE)+'\n'

    input_NMI = ski.metrics.normalized_mutual_information(act_output_arr, input_arr, bins=100)
    input_str+= 'NMI: '+'{:0.5f}'.format(input_NMI)+'\n'
    pred_NMI = ski.metrics.normalized_mutual_information(act_output_arr, pred_output_arr, bins=100)
    pred_str+= 'NMI: '+'{:0.5f}'.format(pred_NMI)+'\n'
    actual_NMI = ski.metrics.normalized_mutual_information(act_output_arr, act_output_arr, bins=100)
    act_str+= 'NMI: '+'{:0.5f}'.format(actual_NMI)+'\n'

    inversed_input_arr = img_inverse_transform(input_arr)
    inversed_pred_output_arr = img_inverse_transform(pred_output_arr)
    inversed_act_output_arr = img_inverse_transform(act_output_arr)
    input_PSNR = ski.metrics.peak_signal_noise_ratio(inversed_act_output_arr, inversed_input_arr, data_range=np.max(inversed_input_arr) - np.min(inversed_input_arr))
    input_str+= 'PSNR: '+'{:0.5f}'.format(input_PSNR)+'\n'
    pred_PSNR = ski.metrics.peak_signal_noise_ratio(inversed_act_output_arr, inversed_pred_output_arr, data_range=np.max(inversed_pred_output_arr) - np.min(inversed_pred_output_arr))
    pred_str+= 'PSNR: '+'{:0.5f}'.format(pred_PSNR)+'\n'
    actual_PSNR = ski.metrics.peak_signal_noise_ratio(inversed_act_output_arr, inversed_act_output_arr, data_range=np.max(inversed_act_output_arr) - np.min(inversed_act_output_arr))
    act_str+= 'PSNR: '+'{:0.5f}'.format(actual_PSNR)+'\n'

    input_SSI = ski.metrics.structural_similarity(act_output_arr, input_arr, gradient=False, data_range=np.max(input_arr) - np.min(input_arr), channel_axis=None)
    input_str+= 'SSI: '+'{:0.5f}'.format(input_SSI)+'\n'
    pred_SSI = ski.metrics.structural_similarity(act_output_arr, pred_output_arr, gradient=False, data_range=np.max(pred_output_arr) - np.min(pred_output_arr), channel_axis=None)
    pred_str+= 'SSI: '+'{:0.5f}'.format(pred_SSI)+'\n'
    actual_SSI = ski.metrics.structural_similarity(act_output_arr, act_output_arr, gradient=False, data_range=np.max(act_output_arr) - np.min(act_output_arr), channel_axis=None)
    act_str+= 'SSI: '+'{:0.5f}'.format(actual_SSI)+'\n'

    input_total_density = np.sum(inversed_input_arr)
    input_str+= 'Total Density: '+'{:0.5f}'.format(input_total_density)+'\n'
    pred_total_density = np.sum(inversed_pred_output_arr)
    pred_str+= 'Total Density: '+'{:0.5f}'.format(pred_total_density)+'\n'
    act_total_density = np.sum(inversed_act_output_arr)
    act_str+= 'Total Density: '+'{:0.5f}'.format(act_total_density)+'\n'

    input_FFTL = F.mse_loss(torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(input_arr)))), torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(act_output_arr)))))
    input_str+= 'FFTL: '+'{:0.5f}'.format(input_FFTL)+'\n'
    pred_FFTL = F.mse_loss(torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(pred_output_arr)))), torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(act_output_arr)))))
    pred_str+= 'FFTL: '+'{:0.5f}'.format(pred_FFTL)+'\n'
    actual_FFTL = F.mse_loss(torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(act_output_arr)))), torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(act_output_arr)))))
    act_str+= 'FFTL: '+'{:0.5f}'.format(actual_FFTL)+'\n'

    fig = plt.figure(figsize=(12, 12))
    
    for i in range(3):
        ax  = fig.add_subplot(3, 3, 3*i+1)
        c = ax.pcolormesh(np.mean(input_arr, axis=i).T, vmin=-1, vmax=1, cmap='gnuplot2')
        fig.colorbar(c, ax=ax)
        if i==0:
            ax.set_title('Input: '+str(label[0])+'_'+str(label[1])+'_'+str(label[2]))
        ax.grid('both')
        if i==2:
            ax.text(0.15, -0.02, input_str, fontsize=10, transform=plt.gcf().transFigure, bbox=dict(boxstyle="round", edgecolor='black', facecolor='white'))

    for i in range(3):
        ax2  = fig.add_subplot(3, 3, 3*i+2)
        c = ax2.pcolormesh(np.mean(pred_output_arr, axis=i).T, vmin=-1, vmax=1, cmap='gnuplot2')
        fig.colorbar(c, ax=ax2)
        if i==0:
            ax2.set_title('Predicted Output '+str(label[0])+'_'+str(label[1])+'_'+str(label[2]))
        ax2.grid('both')
        if i==2:
            ax2.text(0.4, -0.02, pred_str, fontsize=10, transform=plt.gcf().transFigure, bbox=dict(boxstyle="round", edgecolor='black', facecolor='white'))

    for i in range(3):
        ax3  = fig.add_subplot(3, 3, 3*i+3)
        c = ax3.pcolormesh(np.mean(act_output_arr, axis=i).T, vmin=-1, vmax=1, cmap='gnuplot2')
        fig.colorbar(c, ax=ax3)
        if i==0:
            ax3.set_title('Actual Output '+str(label[0])+'_'+str(label[1])+'_'+str(label[2]))
        ax3.grid('both')
        if i==2:
            ax3.text(0.7, -0.02, act_str, fontsize=10, transform=plt.gcf().transFigure, bbox=dict(boxstyle="round", edgecolor='black', facecolor='white'))


    if save_filename:
        plt.savefig(save_filename, bbox_inches='tight')
    
    if dont_show==False:
        plt.show()
    plt.close()


def compare_output_batch(dataset, args, argsGRU, batch_size=1, num_workers=0):
    """
    Compare model predictions with ground truth using batch processing.
    More efficient for large datasets.
    
    Parameters:
    -----------
    dataset : torch.utils.data.Dataset
        Dataset containing input/output pairs
    args : object
        Arguments object containing model configuration
    argsGRU : object
        GRU-specific arguments
    batch_size : int
        Batch size for processing (default: 1 for individual image saving)
    num_workers : int
        Number of DataLoader workers
    
    Returns:
    --------
    tuple : (mean_MSE, mean_density_difference)
    """
    from torch.utils.data import DataLoader
    
    # Create DataLoader
    loader = DataLoader(
        dataset, 
        batch_size=batch_size, 
        shuffle=False,  # Keep original order
        num_workers=num_workers,
        pin_memory=True
    )
    
    # Load the trained model
    model = hybridConvGRU3DNET(args, argsGRU).to(DEVICE, non_blocking=True)
    checkpoint = torch.load(f"./models/{args.run_name}_{MODELFILE}", map_location=DEVICE, weights_only=False)
    model.load_state_dict(checkpoint['model_state_dict'])
    print(f"Model loaded with {count_parameters(model):,} trainable parameters.", flush=True)
    
    # Set the model to evaluation mode
    model.eval()

    pred_MSEs = []
    density_differences = []
    
    with torch.no_grad():
        for batch_idx, batch_data in enumerate(loader):
            # Unpack batch data
            input_batch, output_batch, label_batch = batch_data
            
            # Format for model input
            inputs = input_batch.unsqueeze(2)  # Add channel dimension
            targets = output_batch.unsqueeze(1)  # Add channel dimension
            
            # Move to device
            inputs, targets = inputs.to(DEVICE, non_blocking=True), targets.to(DEVICE, non_blocking=True)
            
            # Get model predictions
            outputs = model(inputs)
            
            # Process each item in the batch
            for i in range(inputs.size(0)):
                # Extract single sample
                single_input = inputs[i:i+1]
                single_output = outputs[i:i+1]
                single_target = targets[i:i+1]
                single_label = [label_batch[0][i].item(), label_batch[1][i], label_batch[2][i].item()]
                
                # Convert to numpy
                output_np = single_output.cpu().numpy()
                target_np = single_target.cpu().numpy()
                input_np = single_input.cpu().numpy()
                
                # Calculate performance metrics
                pred_MSE, density_difference = model_performance(output_np, target_np)
                
                # Create comparison plot
                sample_idx = batch_idx * batch_size + i
                plot_comparison(
                    input_np[0][-1][0],  # Last frame of input sequence
                    output_np[0][0], 
                    target_np[0][0], 
                    single_label, 
                    save_filename=f'./models/plots/output/img_' + "_".join(str(x) for x in single_label) + '.png', 
                    dont_show=True
                )
                
                pred_MSEs.append(pred_MSE)
                density_differences.append(density_difference)
    
    return np.mean(pred_MSEs), np.mean(density_differences)

import time

start = time.time()
print('Preparing HDF5 datasets!', flush=True)
train_dataset, test_dataset, train_indices, test_indices = load_hdf5_data_for_training(DATASETNAME, num_samples=2000,test_percentage=TEST_PERCENTAGE,seed=seed)
end = time.time()
print(f"Time taken to load image hdf5: {end - start:.2f} seconds", flush=True)
start = time.time()
#launch training
parser = argparse.ArgumentParser()
argsUNET = parser.parse_args(args=[])
argsUNET.image_size = IMAGESIZE
argsUNET.device = DEVICE
argsUNET.run_name = "convGRUnet3D_best_GPU_optimized_accumulated"
argsUNET.loss_type = 'dynamic'
argsUNET.epochs = 100
argsUNET.lr = 3e-4
argsUNET.batch_size = 4
argsUNET.accum_steps= 2
argsUNET.num_workers = 4
argsUNET.Adamw_weight_decay = 1e-4
argsUNET.DWA_temperature = 2
argsUNET.base_f = 32
argsUNET.depth = 3
argsUNET.attention_int_division = 2
argsUNET.conv_param = (3,1,1,1)  # kernel size, stride, padding, dilation
argsUNET.LRLU_slope = 0.01  # LeakyReLU slope, 0 for ReLU
argsUNET.dropout_rate = 0.2  # dropout rate, 0 for no dropout
argsGRU = parser.parse_args(args=[])
argsGRU.in_channel = 1
argsGRU.out_channel = 1
argsGRU.in_seq = 2
argsGRU.out_seq = 1
argsGRU.num_layers = 2
argsGRU.kernels = [3,3]
argsGRU.hidden_dims = [32,64]
argsGRU.bias = False

setup_memory_optimization()
print(argsUNET.run_name + " Starting training!", flush=True)
train_unet_with_dataset_optimized_w_accumulation(train_dataset, argsUNET, argsGRU) 
#train_unet(input_image_train, output_image_train, label_train, argsUNET, argsGRU)
print(argsUNET.run_name + " Training completed!", flush=True)
end = time.time()
print(f"Time taken for training {argsUNET.run_name}: {end - start:.2f} seconds", flush=True)
#launch testing
avg_test_MSE, avg_test_density_diff = compare_output_batch(test_dataset, argsUNET, argsGRU, batch_size=2, num_workers=4)

print(f"Model: {argsUNET.run_name} Avg Test MSE: {avg_test_MSE:.4f}, Avg Test Density Difference: {avg_test_density_diff:.4f}", flush=True)
print(f"Model: {argsUNET.run_name} Testing completed!", flush=True)

