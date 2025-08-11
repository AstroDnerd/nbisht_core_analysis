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

DATAFILE  = '/home/x-nbisht1/scratch/projects/128/B2/datasets/nb101_ML_dataset_AllData_AutoEnc.h5'
IMAGEINPUT = '/home/x-nbisht1/scratch/projects/128/B2/datasets/Unet/Images_seq2_Input_3D/'
IMAGEOUTPUT = '/home/x-nbisht1/scratch/projects/128/B2/datasets/Unet/Images_seq2_Output_3D/'
MODELFILE = 'nnmodel.dict'

IMAGESIZE = 64
FRAME_DIFF = 30
TEST_PERCENTAGE = 0.2

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

def compare_output(input_arr, output_arr, labels, args, argsGRU):
    # Load the trained model
    model = hybridConvGRU3DNET(args, argsGRU).to(DEVICE, non_blocking=True)
    model.load_state_dict(torch.load(f"./models/{args.run_name}_{MODELFILE}", map_location=DEVICE))
    print(f"Model loaded with {count_parameters(model):,} trainable parameters.", flush=True)
    # Set the model to evaluation mode
    model.eval()

    pred_MSEs = []
    density_differences = []
    with torch.no_grad():
        for i in range(len(input_arr)):
            inputs = input_arr[i].unsqueeze(0)
            targets = output_arr[i]
            label = labels[i]
            inputs, targets = inputs.to(DEVICE, non_blocking=True), targets.to(DEVICE, non_blocking=True)
            outputs = model(inputs)
            # Save the output images
            output_image = outputs.cpu().numpy()
            pred_MSE, density_difference = model_performance(output_image[0][0], targets[0].numpy())
            pred_MSEs.append(pred_MSE)
            density_differences.append(density_difference)
    return np.mean(pred_MSEs), np.mean(density_differences)

import time

start = time.time()
print('Loading images!', flush=True)
import random
rng_indices = random.sample(range(0, len(os.listdir(IMAGEINPUT))), 2000)  # Randomly select 400 images for training and testing
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
    
end = time.time()
print(f"Total set size: {len(input_image_arr)}, Training set size: {len(input_image_train)}, Test set size: {len(input_image_test)}", flush=True)
print(f"Time taken to load images: {end - start:.2f} seconds", flush=True)
start = time.time()
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
argsGRU.num_layers = 1
argsGRU.kernels = [3]
argsGRU.hidden_dims = [16]
argsGRU.bias = True

print(argsUNET.run_name + " Starting training!", flush=True)
train_unet(input_image_train, output_image_train, label_train, argsUNET, argsGRU)
print(argsUNET.run_name + " Training completed!", flush=True)
end = time.time()
print(f"Time taken for training {argsUNET.run_name}: {end - start:.2f} seconds", flush=True)
#launch testing
avg_test_MSE, avg_test_density_diff = compare_output(input_image_test, output_image_test, label_test, argsUNET, argsGRU)

print(f"Model: {argsUNET.run_name} Testing completed!", flush=True)

