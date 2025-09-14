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
from model5_modules import *

from tqdm import tqdm

# set a seed to ensure reproducibility
seed = 128
rnd  = np.random.RandomState(seed)

DATAFILE  = '/home/x-nbisht1/scratch/projects/128/B2/datasets/nb101_ML_dataset_AllData_AutoEnc.h5'
IMAGEINPUT = '/home/x-nbisht1/scratch/projects/128/B2/datasets/Unet/Images_Input_3D/'
IMAGEOUTPUT = '/home/x-nbisht1/scratch/projects/128/B2/datasets/Unet/Images_Output_3D/'
CORESET  = '/home/x-nbisht1/scratch/projects/128/B2/datasets/nb101_all_frames.h5'
DENSITYCUBES = '/home/x-nbisht1/scratch/projects/128/B2/datasets/nb101_TimeseriesCubes_Density.npy'
MODELFILE = 'nnmodel.dict'

IMAGESIZE = 64
FRAMES = np.arange(0,90, 1)
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

def compare_output(input_arr, output_arr, labels, args):
    # Load the trained model
    model = UNet3D(in_ch=1, out_ch=1, args=args).to(DEVICE)
    checkpoint = torch.load(f"./{args.phase}/{args.run_name}_{MODELFILE}", map_location=DEVICE)
    model.load_state_dict(checkpoint['model_state_dict'])
    print(f"Model loaded for testing with {count_parameters(model):,} trainable parameters.", flush=True)
    # Set the model to evaluation mode
    model.eval()

    pred_MSEs = []
    density_differences = []
    with torch.no_grad():
        for i in range(len(input_arr)):
            inputs = input_arr[i].unsqueeze(0).unsqueeze(0)
            targets = output_arr[i].unsqueeze(0).unsqueeze(0)
            label = labels[i]
            inputs, targets = inputs.to(DEVICE), targets.to(DEVICE)
            outputs = model(inputs)
            # Save the output images
            output_image = outputs.cpu().numpy()+ inputs.cpu().numpy()
            pred_MSE, density_difference = model_performance(output_image[0][0], targets[0][0].numpy())
            pred_MSEs.append(pred_MSE)
            density_differences.append(density_difference)
    return np.mean(pred_MSEs), np.mean(density_differences)

import time

start = time.time()
print('Loading images!', flush=True)
import random
rng_indices = random.sample(range(0, len(os.listdir(IMAGEINPUT))), 400)  # Randomly select 400 images for training and testing
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

model_df = pd.read_csv("./phase1_model_performance.csv")
best_model = model_df.loc[0]

print(best_model, flush=True)
epochs_list = [4,6,8]
for i in range(3):
    start = time.time()
    #launch training
    parser = argparse.ArgumentParser()
    args = parser.parse_args(args=[])
    args.image_size = IMAGESIZE
    args.device = DEVICE
    args.run_name = "Unet_Delta_"+best_model['model_name']
    args.loss_type = 'dynamic'
    args.epochs = epochs_list[i]
    args.lr = float(best_model['learning_rate'])
    args.batch_size = int(best_model['batch_size'])
    args.Adamw_weight_decay = float(best_model['Adamw_weight_decay'])
    args.DWA_temperature = int(best_model['DWA_temperature'])
    args.base_f = int(best_model['base_f'])
    args.depth = int(best_model['depth'])
    args.attention_int_division = int(best_model['attention_int_division'])
    import ast
    args.conv_param = ast.literal_eval(best_model['convolution_param'])  # kernel size, stride, padding, dilation
    args.LRLU_slope = float(best_model['LRLU_slope'])  # LeakyReLU slope, 0 for ReLU
    args.dropout_rate = float(best_model['dropout_rate'])  # dropout rate, 0 for no dropout
    args.phase = 'models_phase2'

    min_loss = train_unet(input_image_train, output_image_train, label_train, args)
    print(best_model['model_name'] + " Training completed!", flush=True)
    end = time.time()
    print(f"Time taken for training {best_model['model_name']}: {end - start:.2f} seconds", flush=True)
    #launch testing
    avg_test_MSE, avg_test_density_diff = compare_output(input_image_test, output_image_test, label_test, args)
    print(f"Average Test MSE: {avg_test_MSE}, Average Test Density Difference: {avg_test_density_diff}", flush=True)
