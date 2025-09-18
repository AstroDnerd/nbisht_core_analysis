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
from model3_0_modules import *

from tqdm import tqdm

# set a seed to ensure reproducibility
seed = 128
rnd  = np.random.RandomState(seed)

IMAGEINPUT = '/home/x-nbisht1/scratch/projects/128/B2/datasets/Unet/Images_seq2_Input_3D/'
IMAGEOUTPUT = '/home/x-nbisht1/scratch/projects/128/B2/datasets/Unet/Images_seq2_Output_3D/'

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


def compare_output(input_arr, output_arr, labels, args, argsGRU):
    # Load the trained model
    model = hybridUNETconvGRU3D(args, argsGRU).to(DEVICE, non_blocking=True)
    checkpoint = torch.load(f"./models/{args.run_name}_{MODELFILE}", map_location=DEVICE, weights_only=False)
    model.load_state_dict(checkpoint['model_state_dict'])
    print(f"Model loaded with {count_parameters(model):,} trainable parameters.", flush=True)
    # Set the model to evaluation mode
    model.eval()

    output_dir = './models/plots/test/output_'+argsUNET.run_name
    make_dir(output_dir)
    pred_MSEs = []
    density_differences = []
    with torch.no_grad():
        for i in range(len(input_arr)):
            inputs = input_arr[i].unsqueeze(0).unsqueeze(2)
            targets = output_arr[i].unsqueeze(0).unsqueeze(0)
            label = labels[i]
            inputs, targets = inputs.to(DEVICE, non_blocking=True), targets.to(DEVICE, non_blocking=True)
            outputs = model(inputs)
            # Save the output images
            output_image = outputs.cpu().numpy()
            pred_MSE, density_difference = model_performance(output_image, targets.numpy())
            plot_comparison(inputs[0][-1][0].numpy(), output_image[0][0], targets[0][0].numpy(), label, save_filename=output_dir+'/img'+"_".join(str(x) for x in label)+'.png', dont_show = True)
            pred_MSEs.append(pred_MSE)
            density_differences.append(density_difference)
    return np.mean(pred_MSEs), np.mean(density_differences)




import time

start = time.time()
print('Loading images!', flush=True)
import random
rng_indices = random.sample(range(0, len(os.listdir(IMAGEINPUT))), 2000)  # Randomly select 2000 images for training and testing
input_image_dic = get_subset_images(IMAGEINPUT, rng_indices)
output_image_dic = get_subset_images(IMAGEOUTPUT, rng_indices)
input_image_arr = []
output_image_arr = []
label_arr = []
for key in input_image_dic.keys():
    input_image_arr.append(torch.tensor(input_image_dic[key]['img']))
    output_image_arr.append(torch.tensor(output_image_dic[key]['img']))
    label_arr.append(input_image_dic[key]['label'])
input_image_train, input_image_test, output_image_train, output_image_test, label_train, label_test = train_test_split(input_image_arr, output_image_arr, label_arr, 
                                                                                                                       shuffle = False, test_size=TEST_PERCENTAGE, random_state=seed)

#validation set:
val_indices = np.delete(range(0, len(os.listdir(IMAGEINPUT))), rng_indices).tolist()
input_image_val_dic = get_subset_images(IMAGEINPUT, val_indices)
output_image_val_dic = get_subset_images(IMAGEOUTPUT, val_indices)
input_val_arr = []
output_val_arr = []
label_val_arr = []
for key in input_image_val_dic.keys():
    input_val_arr.append(torch.tensor(input_image_val_dic[key]['img']))
    output_val_arr.append(torch.tensor(output_image_val_dic[key]['img']))
    label_val_arr.append(input_image_val_dic[key]['label'])


end = time.time()
print(f"Total set size: {len(input_image_arr)}, Training set size: {len(input_image_train)}, Test set size: {len(input_image_test)}, validation set: {len(input_val_arr)}, ", flush=True)
print(f"Time taken to load images: {end - start:.2f} seconds", flush=True)
start = time.time()
#launch training
parser = argparse.ArgumentParser()
argsUNET = parser.parse_args(args=[])
argsUNET.image_size = IMAGESIZE
argsUNET.device = DEVICE
argsUNET.run_name = 'UNetconvGRU3D_TrainingPhases_yes_validation'
argsUNET.loss_type = 'dynamic'
argsUNET.in_channel = 1
argsUNET.out_channel = 1
argsUNET.epochs = 100
argsUNET.lr = 1e-3 #starting lr
argsUNET.batch_size = 8 #best
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

print(argsUNET.run_name + " Starting training!", flush=True)
#train_unet(input_image_train, output_image_train, label_train, argsUNET, argsGRU, validation_input = input_val_arr, validation_output = output_val_arr)
train_unet_training_phase(input_image_train, output_image_train, label_train, argsUNET, argsGRU, validation_input = input_val_arr, validation_output = output_val_arr)
print(argsUNET.run_name + " Training completed!", flush=True)
end = time.time()
print(f"Time taken for training {argsUNET.run_name}: {end - start:.2f} seconds", flush=True)
#launch testing
avg_test_MSE, avg_test_density_diff = compare_output(input_image_test, output_image_test, label_test, argsUNET, argsGRU)

print(f"Model: {argsUNET.run_name} Avg Test MSE: {avg_test_MSE:.4f}, Avg Test Density Difference: {avg_test_density_diff:.4f}", flush=True)
print(f"Model: {argsUNET.run_name} Testing completed!", flush=True)


IMAGEINPUT_u502 = '/home/x-nbisht1/scratch/projects/128/B2/datasets/Unet/u502_Images_seq2_Input_3D/'
IMAGEOUTPUT_u502 = '/home/x-nbisht1/scratch/projects/128/B2/datasets/Unet/u502_Images_seq2_Output_3D/'

print('Loading u502 images!', flush=True)
import random
rng_indices = random.sample(range(0, len(os.listdir(IMAGEINPUT_u502))), 1000)  # Randomly select 2000 images for testing
input_image_dic = get_subset_images(IMAGEINPUT_u502, rng_indices)
output_image_dic = get_subset_images(IMAGEOUTPUT_u502, rng_indices)
input_image_arr = []
output_image_arr = []
label_arr = []
for key in input_image_dic.keys():
    input_image_arr.append(torch.tensor(input_image_dic[key]['img']))
    output_image_arr.append(torch.tensor(output_image_dic[key]['img']))
    label_arr.append(input_image_dic[key]['label'])

def plot_comparison_u502(input_arr, pred_output_arr, act_output_arr, label, loss_dict, save_filename='xy_data.png', dont_show = False):
    input_str = 'INPUT IMAGE:\n'
    pred_str = 'PREDICTED IMAGE:\n'
    act_str = 'ACTUAL IMAGE:\n'
    # Calculate metrics

    input_L1 = F.l1_loss(torch.from_numpy(act_output_arr), torch.from_numpy(input_arr))
    loss_dict['l1']['input'].append(input_L1)
    input_str+= 'L1: '+'{:0.5f}'.format(input_L1)+'\n'
    pred_L1 = F.l1_loss(torch.from_numpy(act_output_arr), torch.from_numpy(pred_output_arr))
    loss_dict['l1']['pred'].append(pred_L1)
    pred_str+= 'L1: '+'{:0.5f}'.format(pred_L1)+'\n'
    actual_L1 = F.l1_loss(torch.from_numpy(act_output_arr), torch.from_numpy(act_output_arr))
    act_str+= 'L1: '+'{:0.5f}'.format(actual_L1)+'\n'

    input_hist_loss = hist_loss(torch.from_numpy(act_output_arr), torch.from_numpy(input_arr), bins=32)
    loss_dict['hist']['input'].append(input_hist_loss)
    input_str+= 'HistLoss: '+'{:0.5f}'.format(input_hist_loss)+'\n'
    pred_hist_loss = hist_loss(torch.from_numpy(act_output_arr), torch.from_numpy(pred_output_arr), bins=32)
    loss_dict['hist']['pred'].append(pred_hist_loss)
    pred_str+= 'HistLoss: '+'{:0.5f}'.format(pred_hist_loss)+'\n'
    actual_hist_loss = hist_loss(torch.from_numpy(act_output_arr), torch.from_numpy(act_output_arr), bins=32)
    act_str+= 'HistLoss: '+'{:0.5f}'.format(actual_hist_loss)+'\n'

    inversed_input_arr = img_inverse_transform(input_arr)
    inversed_pred_output_arr = img_inverse_transform(pred_output_arr)
    inversed_act_output_arr = img_inverse_transform(act_output_arr)
    input_total_density = np.sum(inversed_input_arr)
    pred_total_density = np.sum(inversed_pred_output_arr)
    act_total_density = np.sum(inversed_act_output_arr)
    input_mass_loss = abs(act_total_density - input_total_density) / (act_total_density+ 1e-8)
    loss_dict['total_density']['input'].append(input_mass_loss)
    input_str+= 'Total Density Loss: '+'{:0.5f}'.format(input_mass_loss)+'\n'

    pred_mass_loss = abs(act_total_density - pred_total_density) / (act_total_density+ 1e-8)
    loss_dict['total_density']['pred'].append(pred_mass_loss)
    pred_str+= 'Total Density Loss: '+'{:0.5f}'.format(pred_mass_loss)+'\n'

    act_mass_loss = abs(act_total_density - act_total_density) / (act_total_density+ 1e-8)
    act_str+= 'Total Density Loss: '+'{:0.5f}'.format(act_mass_loss)+'\n'

    input_FFTL = F.mse_loss(torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(input_arr)))), torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(act_output_arr)))))
    loss_dict['fft']['input'].append(input_FFTL)
    input_str+= 'FFTL: '+'{:0.5f}'.format(input_FFTL)+'\n'
    pred_FFTL = F.mse_loss(torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(pred_output_arr)))), torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(act_output_arr)))))
    loss_dict['fft']['pred'].append(pred_FFTL)
    pred_str+= 'FFTL: '+'{:0.5f}'.format(pred_FFTL)+'\n'
    actual_FFTL = F.mse_loss(torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(act_output_arr)))), torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(act_output_arr)))))
    act_str+= 'FFTL: '+'{:0.5f}'.format(actual_FFTL)+'\n'

    hd_input = torch.heaviside(torch.from_numpy(input_arr)-torch.quantile(torch.from_numpy(input_arr.flatten()), 0.99), torch.tensor([1.0]))
    hd_pred = torch.heaviside(torch.from_numpy(pred_output_arr)-torch.quantile(torch.from_numpy(pred_output_arr.flatten()), 0.99), torch.tensor([1.0]))
    hd_output = torch.heaviside(torch.from_numpy(act_output_arr)-torch.quantile(torch.from_numpy(act_output_arr.flatten()), 0.99), torch.tensor([1.0]))
    hd_input_loss = F.l1_loss(hd_input, hd_output).sum().item()
    loss_dict['hd']['input'].append(hd_input_loss)
    input_str+= 'HD Loss: '+'{:0.5f}'.format(hd_input_loss)+'\n'
    hd_pred_loss = F.l1_loss(hd_pred, hd_output).sum().item()
    loss_dict['hd']['pred'].append(hd_pred_loss)
    pred_str+= 'HD Loss: '+'{:0.5f}'.format(hd_pred_loss)+'\n'
    actual_HD =F.l1_loss(hd_output, hd_output).sum().item()
    act_str+= 'HD Loss: '+'{:0.5f}'.format(actual_HD)+'\n'

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

def hist_loss(pred, true, bins=16):
    # flatten
    p = pred.view(-1)
    t =  true.view(-1)
    # compute hist counts (differentiable via torch.histc)
    p_h = torch.histc(p, bins=bins, min=-1, max=1)
    t_h = torch.histc(t, bins=bins, min=-1, max=1)
    p_h /= p_h.sum()
    t_h /= t_h.sum()
    return F.l1_loss(p_h, t_h)


from decimal import Decimal
def compare_output_u502(input_arr, output_arr, labels, args, argsGRU):
    # Load the trained model
    model = hybridUNETconvGRU3D(args, argsGRU).to(DEVICE, non_blocking=True)
    checkpoint = torch.load(f"./models/{args.run_name}_{MODELFILE}", map_location=DEVICE, weights_only=False)
    model.load_state_dict(checkpoint['model_state_dict'])
    print(f"Model loaded with {count_parameters(model):,} trainable parameters.", flush=True)
    # Set the model to evaluation mode
    model.eval()

    base_dir_name = './models/plots/test_u502/output_'+args.run_name
    all_pred_dir_name = base_dir_name+'/all_pred'
    make_dir(all_pred_dir_name)

    L1_loss_dict = {'input': [], 'pred': []}
    hist_loss_dict = {'input': [], 'pred': []}
    total_density_loss_dict = {'input': [], 'pred': []}
    fft_loss_dict = {'input': [], 'pred': []}
    hd_loss_dict = {'input': [], 'pred': []}
    Loss_dict = {'l1': L1_loss_dict, 'hist': hist_loss_dict, 'total_density': total_density_loss_dict, 'fft': fft_loss_dict, 'hd': hd_loss_dict}

    with torch.no_grad():
        for i in range(len(input_arr)):
            inputs = input_arr[i].unsqueeze(0).unsqueeze(2)
            targets = output_arr[i].unsqueeze(0).unsqueeze(0)
            label = labels[i]
            inputs, targets = inputs.to(DEVICE, non_blocking=True), targets.to(DEVICE, non_blocking=True)
            outputs = model(inputs)
            # Save the output images
            output_image = outputs.cpu().numpy()
            plot_comparison_u502(inputs[0][-1][0].numpy(), output_image[0][0], targets[0][0].numpy(), label, Loss_dict, save_filename=all_pred_dir_name+'/img'+"_".join(str(x) for x in label)+'.png', dont_show = True)
    print("Testing completed and output images plotted.", flush=True)

    pred_data = []
    input_data = []
    Loss_names = Loss_dict.keys()
    for key in Loss_names:
        pred_data.append(np.sort(Loss_dict[key]['pred']))
        input_data.append(np.sort(Loss_dict[key]['input']))
    
    pred_data = np.array(pred_data).T
    input_data = np.array(input_data).T

    fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(24, 16))
    parts = ax1.violinplot(
            pred_data, showmeans=False, showmedians=False,
            showextrema=False)

    for pc in parts['bodies']:
        pc.set_facecolor('plum')
        pc.set_edgecolor('indigo')
        pc.set_alpha(1)

    quartile1, medians, quartile3 = np.percentile(pred_data, [25, 50, 75], axis=0)
    medians_orig= np.percentile(input_data,  50, axis=0)
    whiskers = np.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(pred_data, quartile1, quartile3)])
    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]

    inds = np.arange(1, len(medians) + 1)
    #prediction norm
    ax1.scatter(inds, medians, marker='o', color='white', s=30, zorder=3)
    for i, txt in enumerate(medians):
        ax1.annotate('%.3E'%Decimal(txt), (inds[i]+0.1, medians[i]))
    ax1.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
    ax1.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

    #original as predicted
    ax1.scatter(inds, medians_orig, marker='x', color='red', s=50, zorder=3)
    for i, txt in enumerate(medians_orig):
        ax1.annotate('%.3E'%Decimal(txt), (inds[i]-0.1, medians_orig[i]))
    print(medians_orig)

    ax1.axhline(0, color='black', linestyle='--', lw=2)
    ax1.set_ylabel('Loss value', fontsize=18)
    # set style for the axes
    for ax in [ax1]:
        set_axis_style(ax, Loss_names)
    
    ax1.set_xlabel('Losses', fontsize=18)
    ax1.set_yscale('log')
    ax1.set_ylim([1e-10,2])
    ax1.set_title(args.run_name)


    plt.savefig(base_dir_name+'/TestSummary.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


#launch testing
compare_output_u502(input_image_arr, output_image_arr, label_arr, argsUNET, argsGRU)

print(f"Model: {argsUNET.run_name} Diffsim Testing completed!", flush=True)

