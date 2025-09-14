# standard system modules
import os, sys
import h5py 
import argparse
# standard module for tabular data
import json

# standard module for array manipulation
import numpy as np
from itertools import permutations

# standard statistical module
import scipy.stats as st
from scipy import linalg
from scipy.stats import ks_2samp
import skimage as ski
from skimage.metrics import structural_similarity as ski_ssim, normalized_mutual_information as ski_nmi

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
from torch.utils.data import DataLoader, TensorDataset
from torchvision.utils import save_image
from tqdm import tqdm


from model1_0_modules import *

# set a seed to ensure reproducibility
seed = 128
DATAFILE  = '/home/x-nbisht1/scratch/projects/128/B2/datasets/nb101_ML_dataset_AllData_AutoEnc.h5'
IMAGEINPUT = '/home/x-nbisht1/scratch/projects/128/B2/datasets/Unet/Images_seq2_Input_3D/'
IMAGEOUTPUT = '/home/x-nbisht1/scratch/projects/128/B2/datasets/Unet/Images_seq2_Output_3D/'
MODELFILE = 'nnmodel_lr3e4.dict'

IMAGESIZE = 64
FRAME_DIFF = 30
TEST_PERCENTAGE = 0.01

DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
#DEVICE = torch.device('cpu')

print(f'Available device: {str(DEVICE):4s}')


import time

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
    

def plot_comparison(input_arr, pred_output_arr, act_output_arr, label, save_filename='xy_data.png', dont_show = False):
    # Calculate metrics


    input_flat = np.ravel(input_arr)
    pred_flat = np.ravel(pred_output_arr) 
    actual_flat = np.ravel(act_output_arr)


    fig = plt.figure(figsize=(12, 8))
    ax  = fig.add_subplot(3, 1, 1)
    counts, bins, patches = ax.hist(input_flat, bins=16, alpha=0.5, range = (-4, 4))
    for count, bin_edge, patch in zip(counts, bins, patches):
        # Get the x-coordinate (center of the bar) and y-coordinate (top of the bar)
        x = patch.get_x() + patch.get_width() / 2
        y = patch.get_height()

        # Add the text label slightly above the bar
        if count > 0:  # Only display for non-empty bins
            ax.text(x, y + 0.05 * y, int(count), ha='center', va='bottom', fontsize=8)

    ax.set_title('Input: '+str(label[0])+'_'+str(label[1])+'_'+str(label[2]))
    ax.grid('both')


    ax2  = fig.add_subplot(3, 1, 2)
    counts, bins, patches = ax2.hist(pred_flat, bins=16, alpha=0.5, range = (-4, 4))
    ax2.set_title('Predicted Output '+str(label[0])+'_'+str(label[1])+'_'+str(label[2]))
    for count, bin_edge, patch in zip(counts, bins, patches):
        # Get the x-coordinate (center of the bar) and y-coordinate (top of the bar)
        x = patch.get_x() + patch.get_width() / 2
        y = patch.get_height()

        # Add the text label slightly above the bar
        if count > 0:  # Only display for non-empty bins
            ax2.text(x, y + 0.05 * y, int(count), ha='center', va='bottom', fontsize=8)

    ax2.grid('both')

    ax3  = fig.add_subplot(3, 1, 3)
    counts, bins, patches = ax3.hist(actual_flat, bins=16, alpha=0.5, range = (-4, 4))
    for count, bin_edge, patch in zip(counts, bins, patches):
        # Get the x-coordinate (center of the bar) and y-coordinate (top of the bar)
        x = patch.get_x() + patch.get_width() / 2
        y = patch.get_height()

        # Add the text label slightly above the bar
        if count > 0:  # Only display for non-empty bins
            ax3.text(x, y + 0.05 * y, int(count), ha='center', va='bottom', fontsize=8)

    ax3.set_title('Actual Output '+str(label[0])+'_'+str(label[1])+'_'+str(label[2]))
    ax3.grid('both')

    if save_filename:
        plt.savefig(save_filename, bbox_inches='tight')
    
    if dont_show==False:
        plt.show()
    plt.close()


def compare_output(input_arr, output_arr, labels, args, argsGRU):
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
        for i in range(len(input_arr)):
            inputs = input_arr[i].unsqueeze(0).unsqueeze(2)
            targets = output_arr[i].unsqueeze(0).unsqueeze(0)
            label = labels[i]
            inputs, targets = inputs.to(DEVICE, non_blocking=True), targets.to(DEVICE, non_blocking=True)
            outputs = model(inputs)
            # Save the output images
            output_image = outputs.cpu().numpy()
            plot_comparison(inputs[0][-1][0].numpy(), output_image[0][0], targets[0][0].numpy(), label, save_filename='./models/plots/output/img'+"_".join(str(x) for x in label)+'.png', dont_show = True)


#launch training
parser = argparse.ArgumentParser()
argsUNET = parser.parse_args(args=[])
argsUNET.image_size = IMAGESIZE
argsUNET.device = DEVICE
argsUNET.run_name = "convGRUnet3D_best"
argsUNET.loss_type = 'dynamic'
argsUNET.epochs = 30
argsUNET.lr = 3e-4
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

#launch testing
compare_output(input_image_test, output_image_test, label_test, argsUNET, argsGRU)

print(f"Model: {argsUNET.run_name} Testing completed!", flush=True)
