import os, sys
import h5py 
import argparse
import pandas as pd
import numpy as np
import scipy.stats as st
from scipy import linalg
from scipy.stats import ks_2samp
import skimage as ski
from PIL import Image
import matplotlib as mp
import matplotlib.pyplot as plt
mp.rcParams.update(mp.rcParamsDefault)
import imageio.v3 as im
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
from model1_0 import *
from all_modules import *
from tqdm import tqdm
import time

# Import the profiler
from training_profiler import TrainingProfiler, run_pytorch_profiler

# Your existing constants
seed = 128
rnd = np.random.RandomState(seed)
DATASETNAME = '/anvil/scratch/x-nbisht1/projects/512/NonsinkSimSuite/training_copy.hdf5'
MODELFILE = 'nnmodel.dict'
IMAGESIZE = 64
FRAME_DIFF = 33
TEST_PERCENTAGE = 0.01
DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

print(f'Available device: {str(DEVICE):4s}', flush=True)

# Initialize profiler
profiler = TrainingProfiler(output_dir="profiling_results")

# PROFILE DATA LOADING
print('\n' + '='*80)
print('PROFILING: Loading HDF5 datasets')
print('='*80)

with profiler.timer("total_data_loading"):
    with profiler.timer("hdf5_file_open"):
        # Measure just opening the file
        h5file = h5py.File(DATASETNAME, 'r')
        h5file.close()
    
    # Load datasets
    number_of_samples = 1000
    with profiler.timer("dataset_preparation"):
        train_dataset, test_dataset, train_indices, test_indices = load_hdf5_data_for_training(
            DATASETNAME, 
            num_samples=number_of_samples,
            test_percentage=TEST_PERCENTAGE,
            seed=seed
        )

profiler.record_memory("after_data_loading")

# PROFILE MODEL SETUP
print('\n' + '='*80)
print('PROFILING: Model initialization')
print('='*80)

parser = argparse.ArgumentParser()
argsUNET = parser.parse_args(args=[])
argsUNET.image_size = IMAGESIZE
argsUNET.device = 'cuda' if torch.cuda.is_available() else 'cpu'
argsUNET.run_name = "model1_0_profiling"
argsUNET.loss_type = 'statistical'
argsUNET.in_channel = 4
argsUNET.out_channel = 4
argsUNET.epochs = 2
argsUNET.lr = 3e-4
argsUNET.batch_size = 32
argsUNET.accum_steps = 1
argsUNET.num_workers = 32
argsUNET.Adamw_weight_decay = 1e-4
argsUNET.DWA_temperature = 2
argsUNET.base_f = 32
argsUNET.depth = 3
argsUNET.attention_int_division = 2
argsUNET.conv_param = (3,1,1,1)
argsUNET.LRLU_slope = 0.01
argsUNET.dropout_rate = 0.2

argsGRU = parser.parse_args(args=[])
argsGRU.in_channel = 4
argsGRU.out_channel = 4
argsGRU.in_seq = 2
argsGRU.out_seq = 2
argsGRU.num_layers = 2
argsGRU.kernels = [3, 3]
argsGRU.hidden_dims = [32, 64]
argsGRU.bias = False

modelclass = hybridConvGRU3DNET

with profiler.timer("model_instantiation"):
    model = modelclass(argsUNET, argsGRU)

with profiler.timer("model_to_device"):
    model = model.to(DEVICE)

profiler.record_memory("after_model_init")

# Count parameters
total_params = sum(p.numel() for p in model.parameters())
trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
print(f"Total parameters: {total_params:,}")
print(f"Trainable parameters: {trainable_params:,}")

# PROFILE DATALOADER
print('\n' + '='*80)
print('PROFILING: DataLoader creation')
print('='*80)

with profiler.timer("dataloader_setup"):
    from torch.utils.data import DataLoader
    loader = DataLoader(
        train_dataset, 
        batch_size=argsUNET.batch_size,
        shuffle=True,
        num_workers=argsUNET.num_workers,
        pin_memory=(DEVICE.type == 'cuda'),
        persistent_workers=(argsUNET.num_workers > 0)
    )

# PROFILE INDIVIDUAL BATCH OPERATIONS
print('\n' + '='*80)
print('PROFILING: Individual batch operations (first 3 batches)')
print('='*80)

model.train()
optimizer = torch.optim.AdamW(
    model.parameters(), 
    lr=argsUNET.lr, 
    weight_decay=argsUNET.Adamw_weight_decay
)

use_amp = DEVICE.type == 'cuda'
if use_amp:
    scaler = torch.amp.GradScaler()

# Profile first 3 batches in detail
for batch_idx, batch_data in enumerate(loader):
    if batch_idx >= 3:
        break
    
    print(f"\n--- Profiling Batch {batch_idx + 1} ---")
    
    with profiler.timer(f"batch_{batch_idx}_total"):
        # Unpack
        with profiler.timer(f"batch_{batch_idx}_unpack"):
            if len(batch_data) == 3:
                x, y, labels = batch_data
                if argsUNET.in_channel == 1:
                    x = x.unsqueeze(2)
                    y = y.unsqueeze(1)
            else:
                x, y = batch_data
        
        # Transfer to device
        with profiler.timer(f"batch_{batch_idx}_transfer"):
            x = x.to(DEVICE, non_blocking=True)
            y = y.to(DEVICE, non_blocking=True)
        
        profiler.record_memory(f"batch_{batch_idx}_after_transfer")
        
        # Forward pass
        optimizer.zero_grad()
        
        with profiler.timer(f"batch_{batch_idx}_forward"):
            if use_amp:
                with torch.amp.autocast(device_type='cuda'):
                    pred = model(x)
                pred_f = pred.float()
                y_f = y.float()
            else:
                pred = model(x)
                pred_f = pred
                y_f = y
        
        profiler.record_memory(f"batch_{batch_idx}_after_forward")
        
        # Loss computation (using simple MSE for profiling)
        with profiler.timer(f"batch_{batch_idx}_loss"):
            loss = F.mse_loss(pred_f, y_f)
        
        # Backward
        with profiler.timer(f"batch_{batch_idx}_backward"):
            if use_amp:
                scaler.scale(loss).backward()
            else:
                loss.backward()
        
        profiler.record_memory(f"batch_{batch_idx}_after_backward")
        
        # Optimizer step
        with profiler.timer(f"batch_{batch_idx}_optimizer_step"):
            if use_amp:
                scaler.step(optimizer)
                scaler.update()
            else:
                optimizer.step()
        
        profiler.record_memory(f"batch_{batch_idx}_after_step")
        
        if DEVICE.type == 'cuda':
            torch.cuda.synchronize()

# RUN PYTORCH BUILT-IN PROFILER
#run_pytorch_profiler(model, loader, DEVICE, profiler.output_dir, num_batches=5)

# PROFILE FULL TRAINING
# Uncomment this section to profile the full training loop
'''
print('\n' + '='*80)
print('PROFILING: Full training run')
print('='*80)

with profiler.timer("full_training"):
    min_loss = train_unet_unified_profiling(
        modelclass, 
        train_dataset, 
        argsUNET, 
        argsGRU,
    )
    print(f"Training completed! Minimum loss: {min_loss:.4f}", flush=True)


# GENERATE REPORTS
print('\n' + '='*80)
print('GENERATING PROFILING REPORTS')
print('='*80)

profiler.summarize()
profiler.plot_results()
profiler.save_json()

print('\n' + '='*80)
print('PROFILING COMPLETE!')
print('='*80)
print(f"\n Results saved to: {profiler.output_dir}")
print("\n Generated files:")
print("   - profiling_summary.png : Visual breakdown of timing and memory")
print("   - profiling_data.json   : Detailed timing data")
print("   - torch_trace.json      : Chrome trace (view at chrome://tracing)")
print("   - torch_profiler/       : TensorBoard profiler data")
print("\n Nikhil Notes:")
print("   1. Check which operations take the most time in the summary")
print("   2. Look for memory spikes in the plots")
print("   3. Examine the Chrome trace for GPU kernel-level detail")
print("   4. Use TensorBoard for interactive profiling:")
print("      tensorboard --logdir=profiling_results/torch_profiler")
'''