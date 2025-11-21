import os
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
import torch.nn.functional as F
from tqdm import tqdm
from torch import optim
import torchvision
from torch.utils.data import DataLoader, TensorDataset, Subset
import json
from torchvision.utils import save_image
import logging
from torch.utils.tensorboard import SummaryWriter
import time
import datetime
logging.basicConfig(format="%(asctime)s - %(levelname)s: %(message)s", level=logging.INFO, datefmt="%I:%M:%S")

DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
#DEVICE = torch.device('cpu')

MODELFILE = 'nnmodel.dict'

def make_dir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def count_parameters(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)

def img_transform(np_arr):
    return np.log10(np_arr)

def img_inverse_transform(np_arr):
    return np.power(10,np_arr)

def delete_folder_contents(folder):
    import shutil
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))



from torch.utils.data import Dataset, DataLoader
import h5py
# Custom Dataset class for HDF5 data
import torch
from torch.utils.data import Dataset
import h5py
import os

class HDF5Dataset(Dataset):
    """
    Custom PyTorch Dataset for loading data from HDF5 files efficiently.
    This allows for lazy loading - data is only loaded when needed.
    Handles multiprocessing safely by re-opening files in each worker.
    """
    def __init__(self, hdf5_file, indices=None):
        self.hdf5_file = hdf5_file
        self.archive = None 
        self.last_pid = None # Track which process opened the file
        
        # Get available indices from the file (Open/Close immediately)
        with h5py.File(hdf5_file, 'r') as f:
            all_keys = list(f['input_images'].keys())
            
        if indices is not None:
            self.keys = [all_keys[i] for i in indices if i < len(all_keys)]
        else:
            self.keys = all_keys
            
    def __len__(self):
        return len(self.keys)
    
    def __getitem__(self, idx):
        # 1. Check process ID to handle DataLoader workers safely
        curr_pid = os.getpid()
        
        # If we haven't opened the file, OR if we are in a new process (worker)
        if self.archive is None or self.last_pid != curr_pid:
            # Close existing handle if it exists (to be safe)
            if self.archive is not None:
                try:
                    self.archive.close()
                except:
                    pass
            
            # Open file and store the current process ID
            self.archive = h5py.File(self.hdf5_file, 'r')
            self.last_pid = curr_pid

        # 2. Access data DIRECTLY
        f = self.archive
        
        try:
            key = self.keys[idx]
            
            # Load input and output images
            input_img = torch.tensor(f['input_images'][key][...], dtype=torch.float32)
            output_img = torch.tensor(f['output_images'][key][...], dtype=torch.float32)
            
            # Load label
            label_str = f['labels'][key][()]
            if hasattr(label_str, 'decode'):
                label_str = label_str.decode('utf-8')
            label_str = eval(label_str)
            
            return input_img, output_img, label_str
            
        except Exception as e:
            # Fallback: If the file handle went stale for some other reason, try one reset
            print(f"Error loading index {idx}, attempting reset: {e}")
            self.archive.close()
            self.archive = h5py.File(self.hdf5_file, 'r')
            f = self.archive
            
            key = self.keys[idx]
            input_img = torch.tensor(f['input_images'][key][...], dtype=torch.float32)
            output_img = torch.tensor(f['output_images'][key][...], dtype=torch.float32)
            label_str = f['labels'][key][()]
            if hasattr(label_str, 'decode'):
                label_str = label_str.decode('utf-8')
            label_str = eval(label_str)
            return input_img, output_img, label_str

    def __del__(self):
        """Cleanup file handle on destruction"""
        if self.archive is not None:
            try:
                self.archive.close()
            except:
                pass

class ChunkedHDF5Manager:
    def __init__(self, hdf5_path, chunk_size=2000, test_ratio=0.1, seed=128, total_samples=None):
        """
        Manages splitting HDF5 keys into chunks and loading them into RAM.
        """
        self.hdf5_path = hdf5_path
        self.chunk_size = chunk_size
        self.seed = seed
        self.total_samples = total_samples

        print(f"Initializing Chunk Manager for {hdf5_path}...", flush=True)
        
        # 1. Get all keys once
        with h5py.File(hdf5_path, 'r') as f:
            self.all_keys = list(f['input_images'].keys())
        
        # 2. Shuffle and Split
        rng = np.random.RandomState(seed)
        rng.shuffle(self.all_keys)
        if self.total_samples is not None:
            self.all_keys = self.all_keys[:self.total_samples]
        n_total = len(self.all_keys)
        n_test = int(n_total * test_ratio)
        
        self.test_keys = self.all_keys[:n_test]
        self.train_keys = self.all_keys[n_test:]
        
        print(f"Total samples: {n_total}. Train: {len(self.train_keys)}, Test: {len(self.test_keys)}")
        
        # 3. Create chunks for training
        self.train_chunks = [
            self.train_keys[i : i + chunk_size] 
            for i in range(0, len(self.train_keys), chunk_size)
        ]
        print(f"Created {len(self.train_chunks)} training chunks of size ~{chunk_size}.")

    def load_data_to_ram(self, keys, desc="Data"):
        """
        Loads a specific list of keys entirely into CPU RAM.
        Returns a TensorDataset.
        """
        inputs = []
        outputs = []
        labels = []
        
        print(f"Loading {len(keys)} samples into RAM ({desc})...", flush=True)
        
        with h5py.File(self.hdf5_path, 'r') as f:
            # Pre-fetching datasets to avoid lookups in loop
            dset_in = f['input_images']
            dset_out = f['output_images']
            dset_lbl = f['labels']
            
            for key in keys:
                # Read numpy arrays
                inputs.append(dset_in[key][...])
                outputs.append(dset_out[key][...])
                
                # Handle label decoding
                lbl = dset_lbl[key][()]
                if hasattr(lbl, 'decode'):
                    lbl = lbl.decode('utf-8')
                labels.append(eval(lbl))
        
        # Convert to Tensor Stack
        # This is the heavy memory operation
        tensor_x = torch.tensor(np.stack(inputs), dtype=torch.float32)
        tensor_y = torch.tensor(np.stack(outputs), dtype=torch.float32)
        # Assuming labels are simple enough to keep as list or convert if needed
        # If labels are not used in loss, we can keep them as is.
        # Custom dataset wrapper might be needed if labels are complex types
        
        return CustomTensorDataset(tensor_x, tensor_y, labels)

class CustomTensorDataset(torch.utils.data.Dataset):
    """Simple wrapper to handle tensors + list of labels"""
    def __init__(self, x, y, labels):
        self.x = x
        self.y = y
        self.labels = labels
        
    def __len__(self):
        return len(self.x)
        
    def __getitem__(self, idx):
        return self.x[idx], self.y[idx], self.labels[idx]

class EarlyStopping:
    """Stop training when validation loss stops improving."""
    
    def __init__(self, patience=3, min_delta=0.001, mode='min'):
        """
        Args:
            patience: Number of epochs to wait
            min_delta: Minimum change to qualify as improvement
            mode: 'min' or 'max'
        """
        self.patience = patience
        self.min_delta = min_delta
        self.mode = mode
        self.counter = 0
        self.best_loss = None
        self.should_stop = False
    
    def __call__(self, val_loss):
        """
        Check if should stop training.
        
        Args:
            val_loss: Current validation loss
        
        Returns:
            bool: True if should stop
        """
        if self.best_loss is None:
            self.best_loss = val_loss
            return False
        
        if self.mode == 'min':
            improved = val_loss < (self.best_loss - self.min_delta)
        else:
            improved = val_loss > (self.best_loss + self.min_delta)
        
        if improved:
            self.best_loss = val_loss
            self.counter = 0
        else:
            self.counter += 1
            if self.counter >= self.patience:
                print(f"\nEarly stopping triggered after {self.counter} epochs without improvement")
                self.should_stop = True
                return True
        
        return False

import random
from sklearn.model_selection import train_test_split
def load_hdf5_data_for_training(hdf5_file, num_samples=2000, test_percentage=0.2, seed=128):
    """
    Load data from HDF5 file and prepare for training
    
    Parameters:
    -----------
    hdf5_file : str
        Path to HDF5 file
    num_samples : int
        Number of samples to randomly select
    test_percentage : float
        Percentage of data to use for testing
    seed : int
        Random seed for reproducibility
        
    Returns:
    --------
    tuple : (train_dataset, test_dataset, train_indices, test_indices)
    """
    # Set random seed
    random.seed(seed)
    np.random.seed(seed)
    
    # Get total number of samples in HDF5 file
    with h5py.File(hdf5_file, 'r') as f:
        total_samples = len(f['input_images'].keys())
    
    print(f'Total samples in HDF5 file: {total_samples}')
    
    # Randomly select indices
    if num_samples > total_samples:
        print(f'Warning: Requested {num_samples} samples but only {total_samples} available. Using all samples.')
        selected_indices = list(range(total_samples))
    else:
        selected_indices = random.sample(range(total_samples), num_samples)
    
    # Split indices into train and test
    train_indices, test_indices = train_test_split(
        selected_indices, 
        test_size=test_percentage, 
        random_state=seed,
        shuffle=False  # Keep same as your original code
    )
    
    print(f'Training samples: {len(train_indices)}, Test samples: {len(test_indices)}')
    
    # Create datasets
    train_dataset = HDF5Dataset(hdf5_file, train_indices)
    test_dataset = HDF5Dataset(hdf5_file, test_indices)
    
    return train_dataset, test_dataset, train_indices, test_indices


def soft_histogram_1d(x, bins=32, vmin=-1.0, vmax=1.0, sigma=None, sample_frac=1.0, max_samples=int(3e5)):
    """
    Differentiable soft histogram for 1D data.
    x : tensor shape [N] (flattened). Accepts float32/float64.
    bins, vmin, vmax: histogram bins and range.
    sigma: Gaussian width for soft-assignment. If None, uses bin width/2.
    sample_frac: fraction of voxels to sample for speed (random).
    max_samples: upper cap of samples to avoid huge memory.
    Returns normalized histogram tensor [bins].
    """
    if x.numel() == 0:
        return torch.zeros(bins, device=x.device, dtype=torch.float32)

    x = x.view(-1)
    N = x.numel()

    # Subsample if requested
    if sample_frac < 1.0 or N > max_samples:
        keep = min(int(N * sample_frac), max_samples)
        if keep < 1:
            keep = 1
        idx = torch.randperm(N, device=x.device)[:keep]
        x = x[idx]
        N = x.numel()

    bins_edges = torch.linspace(vmin, vmax, steps=bins+1, device=x.device, dtype=torch.float32)
    centers = 0.5 * (bins_edges[:-1] + bins_edges[1:])          # [bins]
    bin_width = (vmax - vmin) / bins
    if sigma is None:
        sigma = bin_width * 0.5
    # x: [N,1], centers: [1,bins]
    x_col = x.view(-1,1).to(torch.float32)
    c_row = centers.view(1,-1)
    diff = (x_col - c_row) / (sigma + 1e-12)    # [N, bins]
    weights = torch.exp(-0.5 * diff*diff)       # gaussian
    hist = weights.sum(dim=0)                   # [bins]
    norm = hist.sum()
    if norm <= 0:
        return torch.zeros_like(hist)
    return hist / (norm + 1e-12)


def spectral_loss_per_channel(pred_c, targ_c, eps=1e-12, sample_frac=1.0, max_samples=None):
    """
    pred_c, targ_c: tensors of shape [B, D, H, W] (single channel)
    returns scalar MSE between log1p magnitudes of FFTs averaged over batch.
    Work in float32 internally; supports subsampling in spatial domain if requested.
    """
    B = pred_c.shape[0]
    loss = 0.0
    for b in range(B):
        p = pred_c[b].to(torch.float32)
        t = targ_c[b].to(torch.float32)
        if max_samples is not None:
            # subsample voxels for speed: reshape to vector, choose subset, put back to approximate FFT?
            # For spectral fidelity better to compute full FFT; use subsampling only if necessary.
            pass
        # full 3D FFT
        p_fft = torch.fft.fftn(p, dim=(-3,-2,-1))
        t_fft = torch.fft.fftn(t, dim=(-3,-2,-1))
        p_pow = torch.log1p(torch.abs(p_fft))
        t_pow = torch.log1p(torch.abs(t_fft))
        loss += F.mse_loss(p_pow, t_pow)
    return loss / float(B)


def soft_hd_mask(x, q=0.999, smooth=1e-2):
    """
    x: [B, ...] tensor
    returns sigmoid((x - thr)/smooth) where thr is the q-quantile per-batch.
    """
    B = x.shape[0]
    masks = []
    flat = x.view(B, -1)
    for b in range(B):
        arr = flat[b]
        thr = torch.quantile(arr, q)
        mask = torch.sigmoid((x[b] - thr) / (smooth + 1e-12))
        masks.append(mask)
    return torch.stack(masks, dim=0)

def mass_and_momentum_losses(pred_density_log, true_density_log, pred_vel=None, true_vel=None, density_is_log=True):
    """
    pred_density_log etc: [B, D, H, W] (single-channel)
    pred_vel, true_vel: [B, 3, D, H, W] or None
    Returns (mass_loss_scalar, momentum_loss_scalar)
    mass computed in linear units assuming density inputs are log10; change if natural log.
    """
    device = pred_density_log.device
    if density_is_log:
        pred_rho = torch.pow(10.0, pred_density_log)
        true_rho = torch.pow(10.0, true_density_log)
    else:
        pred_rho = pred_density_log
        true_rho = true_density_log

    mass_pred = pred_rho.view(pred_rho.shape[0], -1).sum(dim=1)  # [B]
    mass_true = true_rho.view(true_rho.shape[0], -1).sum(dim=1)
    mass_rel = torch.abs(mass_pred - mass_true) / (mass_true + 1e-12)
    mass_loss = torch.log1p(mass_rel).mean()

    momentum_loss = torch.tensor(0.0, device=device)
    if (pred_vel is not None) and (true_vel is not None):
        # pred_vel: [B,3,D,H,W]; compute momentum vector = ∫ rho * v dV
        # Expand rho to multiply
        pr = pred_rho
        tr = true_rho
        # ensure shapes
        # multiply per-component
        mom_pred = (pr.unsqueeze(1) * pred_vel).view(pred_vel.shape[0], pred_vel.shape[1], -1).sum(dim=2)
        mom_true = (tr.unsqueeze(1) * true_vel).view(true_vel.shape[0], true_vel.shape[1], -1).sum(dim=2)
        mom_rel = torch.abs(mom_pred - mom_true) / (torch.abs(mom_true) + 1e-12)
        momentum_loss = mom_rel.mean()
    return mass_loss, momentum_loss

def return_total_loss_multichannel(pred, target, config, input_state=None):
    """
    pred, target: either [B,C,D,H,W] OR [B,T,C,D,H,W]
    config: dict with keys:
      - 'type': prediction type 'normal' or 'delta'
      - 'channel_scales': list of length C (e.g., stds)
      - 'channel_weights': list length C
      - 'density_channel': int (index)
      - 'velocity_channels': list or slice for velocities, e.g. [1,2,3] or None
      - hist params: 'hist_bins','hist_vmin','hist_vmax','hist_sigma','hist_sample_frac'
      - spectral params: 'spec_sample_frac' (optional)
      - hd params: 'hd_q', 'hd_smooth'
      - 'density_is_log' boolean
    Returns: all losses
    """
    # Check prediction type
    pred_type = config.get('type', 'normal')

    # Normalize shape to list of frames
    single_time = False
    if pred.dim() == 5:
        # [B,C,D,H,W] -> add time axis
        pred = pred.unsqueeze(1)
        target = target.unsqueeze(1)
        single_time = True

    # shapes: [B, T, C, D, H, W]
    B, T, C, D, H, W = pred.shape
    device = pred.device

    # Fetch config defaults
    channel_scales = config.get('channel_scales', [1.0]*C)
    channel_weights = config.get('channel_weights', [1.0]*C)
    density_idx = config.get('density_channel', 0)
    vel_idxs = config.get('velocity_channels', None)
    density_is_log = config.get('density_is_log', True)

    # hist params
    hist_bins = config.get('hist_bins', 32)
    hist_vmin = config.get('hist_vmin', -1.0)
    hist_vmax = config.get('hist_vmax', 1.0)
    hist_sigma = config.get('hist_sigma', None)
    hist_sample_frac = config.get('hist_sample_frac', 0.2)

    # spectral params
    spec_sample_frac = config.get('spec_sample_frac', 1.0)
    spec_eval = config.get('do_spectral', True)

    # hd params
    hd_q = config.get('hd_q', 0.99)
    hd_smooth = config.get('hd_smooth', 1e-2)

    # accumulate
    l1_term = torch.tensor(0.0, device=device)
    hist_term = torch.tensor(0.0, device=device)
    spec_term = torch.tensor(0.0, device=device)
    hd_loss_term = torch.tensor(0.0, device=device)
    mass_n_mom_term = torch.tensor(0.0, device=device)

    for t in range(T):
        p_t = pred[:, t]   # [B, C, D, H, W]
        y_t = target[:, t]

        # For delta mode, we need the input state to reconstruct absolute fields
        if pred_type == 'delta':
            if input_state is None:
                raise ValueError("input_state must be provided when config['type']='delta'")
            x_t = input_state[:, t]  # [B, C, D, H, W]

        #Channel-normalized L1
        for c in range(C):
            pc = p_t[:, c].to(torch.float32)
            yc = y_t[:, c].to(torch.float32)
            s = float(channel_scales[c]) if c < len(channel_scales) else 1.0
            wch = float(channel_weights[c]) if c < len(channel_weights) else 1.0
            l1_term = l1_term + wch * F.l1_loss(pc / (s + 1e-12), yc / (s + 1e-12))

        #Histogram loss (per-channel, float32) - do in CPU/float32/autocast-disabled
        with torch.amp.autocast(device_type=str(device), enabled=False):
            for c in range(C):
                pc = p_t[:, c].detach().to(torch.float32).flatten()
                yc = y_t[:, c].detach().to(torch.float32).flatten()
                # use soft histogram; runs on device
                hp = soft_histogram_1d(pc, bins=hist_bins, vmin=hist_vmin, vmax=hist_vmax,
                                       sigma=hist_sigma, sample_frac=hist_sample_frac)
                ht = soft_histogram_1d(yc, bins=hist_bins, vmin=hist_vmin, vmax=hist_vmax,
                                       sigma=hist_sigma, sample_frac=hist_sample_frac)
                hist_term = hist_term + F.l1_loss(hp, ht)
        hist_term = hist_term / float(C)

        #Spectral loss (per-channel) — compute in float32 / autocast disabled
        if spec_eval:
            with torch.amp.autocast(device_type=str(device), enabled=False):
                for c in range(C):
                    pc = p_t[:, c].to(torch.float32)
                    yc = y_t[:, c].to(torch.float32)
                    spec_term = spec_term + spectral_loss_per_channel(pc, yc, sample_frac=spec_sample_frac)
            spec_term = spec_term / float(C)

        #High density (soft mask) — compute only on density channel
        if pred_type == 'delta':
            # Reconstruct absolute density field: density_t+1 = density_t + delta
            # Since density is in log space: log(rho_t+1) = log(rho_t) + delta
            pred_density_reconstructed = x_t[:, density_idx] + p_t[:, density_idx]
            true_density_reconstructed = x_t[:, density_idx] + y_t[:, density_idx]
            
            pd = pred_density_reconstructed.to(torch.float32)
            yd = true_density_reconstructed.to(torch.float32)
        else:
            # Normal mode: use predictions directly
            pd = p_t[:, density_idx].to(torch.float32)
            yd = y_t[:, density_idx].to(torch.float32)
        hd_p = soft_hd_mask(pd, q=hd_q, smooth=hd_smooth)
        hd_y = soft_hd_mask(yd, q=hd_q, smooth=hd_smooth)
        hd_loss_term = hd_loss_term+ F.l1_loss(hd_p, hd_y)

        #Mass & momentum
        vel_present = (vel_idxs is not None) and (len(vel_idxs) >= 3)
        if pred_type == 'delta':
            # Reconstruct density (already done above)
            pred_density_abs = pred_density_reconstructed
            true_density_abs = true_density_reconstructed
            
            # Reconstruct velocity: v_t+1 = v_t + delta_v
            if vel_present:
                pred_vel_reconstructed = x_t[:, vel_idxs] + p_t[:, vel_idxs]
                true_vel_reconstructed = x_t[:, vel_idxs] + y_t[:, vel_idxs]
                
                pred_vel = pred_vel_reconstructed.to(torch.float32)
                true_vel = true_vel_reconstructed.to(torch.float32)
            else:
                pred_vel = None
                true_vel = None
        else:
            # Normal mode: use predictions directly
            pred_density_abs = p_t[:, density_idx]
            true_density_abs = y_t[:, density_idx]
            
            if vel_present:
                pred_vel = p_t[:, vel_idxs].to(torch.float32)
                true_vel = y_t[:, vel_idxs].to(torch.float32)
            else:
                pred_vel = None
                true_vel = None
        mass_l, mom_l = mass_and_momentum_losses(pred_density_abs, true_density_abs, pred_vel, true_vel, density_is_log=density_is_log)
        mass_n_mom_term = mass_n_mom_term + mass_l + mom_l

    # average across time steps
    return l1_term / float(T), hist_term / float(T), spec_term / float(T), hd_loss_term / float(T), mass_n_mom_term / float(T)

def iqr(tensor):
    """
    Calculates the Interquartile Range (IQR) of a PyTorch channel tensor.
    """
    IQR = np.array([])
    for ch in range(tensor.shape[0]):
        # Calculate the 25th percentile (Q1)
        q1 = torch.quantile(tensor[ch], 0.25)
        # Calculate the 75th percentile (Q3)
        q3 = torch.quantile(tensor[ch], 0.75)
        # Calculate the IQR
        IQR = np.append(IQR,float(q3 - q1))
    return IQR

def compute_initial_weights_from_sample(sample_pairs, density_channel=0, hd_q=0.999,
                                        hist_bins=32, hist_sample_frac=0.20, device=None):
    """
    sample_pairs: iterable of tuples (input_tensor, target_tensor)
        where input_tensor has shape [T_in, C, D, H, W] or [C,D,H,W]
        and target_tensor is [C,D,H,W]
    Returns: dict of recommended weights for losses:
         {'w_l1':..., 'w_hist':..., 'w_mass':..., 'w_spec':..., 'w_hd':...}
    Uses persistence predictor: predict = last input frame.
    """

    loss_list = [0.0, 0.0, 0.0, 0.0, 0.0]  # each element: (L1, Hist, Mass, Spec, HD)
    device = sample_pairs[0][0].device if device is None else device
    std_sum = np.array([0.0,0.0,0.0,0.0])
    pair_len = 0
    for (inp, out, labels) in sample_pairs:
        # ensure tensors are on device
        inp = inp.to(device)
        out = out.to(device)
        pair_len+=1
        # choose persistence predictor: last input time frame
        if inp.dim() == 5:
            #[T_in, C, D, H, W]
            y_inp = inp[-1]
            y_pred = out[-1]
        else:
            raise ValueError("unexpected input shape: " + str(inp.shape))

        std_sum+=iqr(y_inp)
        std_sum+=iqr(y_pred)
        std_sum = std_sum/2
        config = {
            'channel_scales': std_sum/pair_len,
            'channel_weights': [1.0, 1.0, 1.0, 1.0],
            'density_channel': density_channel,
            'velocity_channels': [1,2,3],
            'density_is_log': True,
            'hist_bins': hist_bins, 'hist_vmin': -5.0, 'hist_vmax': 5.0, 'hist_sigma': None, 'hist_sample_frac':hist_sample_frac,
            'do_spectral': True, 'spec_sample_frac': 1.0,
            'hd_q': hd_q, 'hd_smooth': 1e-2
        }
        with torch.no_grad():
            l1, hist_l, spectral_l, hd_l, mass_l = return_total_loss_multichannel(inp.unsqueeze(0), out.unsqueeze(0), config)
            loss_list[0]+=l1.item()
            loss_list[1]+=hist_l.item()
            loss_list[2]+=mass_l.item()
            loss_list[3]+=spectral_l.item()
            loss_list[4]+=hd_l.item()
    std_sum = std_sum/pair_len
    return pair_len/np.array(loss_list) , std_sum

def print_memory_stats():
    """Print GPU memory statistics if available"""
    if torch.cuda.is_available():
        print(f"GPU Memory - Allocated: {torch.cuda.memory_allocated()/1024**3:.2f} GB, "
              f"Reserved: {torch.cuda.memory_reserved()/1024**3:.2f} GB", flush=True)

def setup_memory_optimization():
    """
    Set up environment variables for memory optimization
    """
    os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'expandable_segments:True'
    os.environ['CUDA_LAUNCH_BLOCKING'] = '1'  # For debugging
    
    # Additional PyTorch settings
    torch.cuda.empty_cache()
    if torch.cuda.is_available():
        torch.cuda.set_per_process_memory_fraction(0.95)  # Use 95% of available memory
        print(f"GPU Memory: {torch.cuda.get_device_properties(0).total_memory / 1024**3:.1f} GB")
        print(f"GPU Memory allocated: {torch.cuda.memory_allocated() / 1024**3:.1f} GB")
        print(f"GPU Memory reserved: {torch.cuda.memory_reserved() / 1024**3:.1f} GB")

def setup_cpu_optimizations():
    """
    Enable all CPU-specific optimizations.
    
    Args:
        num_threads: Number of threads per operation (4-8 recommended)
    """
    import psutil
    print("SETTING UP CPU OPTIMIZATIONS")
    
    #Enable MKL-DNN (Intel Math Kernel Library)
    torch.backends.mkldnn.enabled = True
    print(f"MKL-DNN enabled: {torch.backends.mkldnn.is_available()}")
    
    # Disable determinism for speed
    torch.backends.mkldnn.deterministic = False
    print(f"Non-deterministic mode (faster)")
    
    #Check bfloat16 support
    use_bf16 = False
    try:
        if hasattr(torch.cpu, 'is_bf16_supported'):
            use_bf16 = torch.cpu.is_bf16_supported()
        elif hasattr(torch.backends, 'mkldnn') and torch.backends.mkldnn.is_available():
            # Try to use bf16
            test = torch.randn(1, 1).bfloat16()
            use_bf16 = True
    except:
        pass
    
    if use_bf16:
        print(f"bfloat16 supported (will use for training)")
    else:
        print(f"bfloat16 not supported (using float32)")
    
    #Enable JIT fusion
    torch.jit.enable_onednn_fusion(True)
    print(f"JIT fusion enabled")
    
    print(f"\nCPU info:")
    print(f"  - Physical cores: {psutil.cpu_count(logical=False)}")
    print(f"  - Logical cores: {psutil.cpu_count(logical=True)}")
    print(f"  - Total RAM: {psutil.virtual_memory().total / 1024**3:.1f} GB")
    
    return use_bf16

def fmt_s(sec):
    return str(datetime.timedelta(seconds=int(sec)))

def print_progress(batch_idx, loader, epoch, start_epoch, args,
                   total_loss, batch_start, training_start,
                   batch_time_ema, ema_alpha=0.15):
    """
    Print an inline status line (overwrites) with ETA, batch timing, GPU memory.
    Returns updated batch_time_ema.

    Required inputs:
      - batch_idx: current batch index (0-based)
      - loader: the DataLoader object
      - epoch: current epoch (0-based or 1-based)
      - start_epoch: epoch resumed from (0 if fresh run)
      - args: args namespace (must contain args.epochs)
      - total_loss: current batch's total_loss (tensor or scalar)
      - batch_start: time.perf_counter() value taken at top of batch loop
      - training_start: time.perf_counter() taken when training run began
      - batch_time_ema: previous EMA value (None to initialize)
      - ema_alpha: smoothing factor (0-1), Larger = faster response.
    """
    # synchronize for accurate timing on GPU
    if torch.cuda.is_available():
        torch.cuda.synchronize()

    batch_time = time.perf_counter() - batch_start
    if batch_time_ema is None:
        batch_time_ema = batch_time
    else:
        batch_time_ema = ema_alpha * batch_time + (1.0 - ema_alpha) * batch_time_ema

    batches_done = batch_idx + 1
    try:
        batches_total = len(loader)
    except Exception:
        batches_total = getattr(loader, "_len_hint", 0) or 0
    batches_left = max(0, batches_total - batches_done)
    #avoid printing every batch if dataloader is very large
    if (batches_done % int(0.1*batches_total)) != 0:
        return batch_time_ema

    est_remain_epoch_s = batch_time_ema * batches_left

    # compute average epoch time so far (approx)
    completed_epochs = max(1, (epoch - start_epoch + 1))
    elapsed_since_training_start = time.perf_counter() - training_start
    avg_epoch_time = elapsed_since_training_start / completed_epochs
    epochs_left = max(0, args.epochs - epoch - 1)
    est_total_remaining_s = est_remain_epoch_s + epochs_left * avg_epoch_time

    # GPU memory stats
    if torch.cuda.is_available():
        mem_alloc_gb = torch.cuda.memory_allocated() / (1024**3)
        mem_max_gb   = torch.cuda.max_memory_allocated() / (1024**3)
        gpu_mem_str = f" | GPU mem {mem_alloc_gb:4.2f}G (peak {mem_max_gb:4.2f}G)"
    else:
        gpu_mem_str = ""

    # current loss
    try:
        current_loss = float(total_loss.item())
    except Exception:
        try:
            current_loss = float(total_loss)
        except Exception:
            current_loss = 0.0

    # Compose inline status
    status = (f"Epoch {epoch+1}/{args.epochs} "
              f"| Batch {batches_done}/{batches_total} "
              f"| loss {current_loss:.4e} "
              f"| batch {batch_time_ema:.3f}s "
              f"| ETA_epoch {fmt_s(est_remain_epoch_s)} "
              f"| ETA_total {fmt_s(est_total_remaining_s)}"
              f"{gpu_mem_str}")

    # print inline (carriage return overwrites)
    print(status + " " * 10, end="\r", flush=True)

    # occasional newline to make log searchable (every 500 batches)
    if (batches_done % 100) == 0:
        print()

    return batch_time_ema


def train_unet_unified_w_mixed_loss(modelclass, dataset, args, argsGRU, validation_dataset=None):
    """
    Optimized training function with mixed precision and improved memory management and gradient accumulation.
    Notes:
    - HDF5 datasets
    - Both CPU and GPU
    - Mixed precision training (GPU only)
    - Gradient accumulation
    - Training phases with adaptive loss scaling
    - Validation
    - Learning rate scheduling
    - Periodic checkpointing
    
    Parameters:
    -----------
    dataset : torch.utils.data.Dataset
        Training dataset (HDF5Dataset)
    args : argparse.Namespace
        Training arguments
    argsGRU : argparse.Namespace
        GRU model arguments
    validation_dataset : torch.utils.data.Dataset, optional
        Validation dataset
    """
    make_dir(f"models/plots/")
    #Setup device
    DEVICE = torch.device(args.device if torch.cuda.is_available() and 'cuda' in args.device else 'cpu')
    print(f"Using device: {DEVICE}", flush=True)
    
    #Memory optimization for GPU
    if DEVICE.type == 'cuda':
        torch.cuda.empty_cache()
        torch.backends.cudnn.benchmark = False
        print_memory_stats()
    
    print(f"Training convGRUNet with {len(dataset)} samples, {args.run_name} run", flush=True)
    
    #DataLoader setup
    num_workers = args.num_workers if hasattr(args, 'num_workers') else min(4, os.cpu_count()//2)
    print(f"Using {num_workers} DataLoader workers", flush=True)
    
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, num_workers=num_workers, pin_memory=(DEVICE.type == 'cuda'), persistent_workers=(num_workers > 0))
    
    # Validation loader if provided
    if validation_dataset is not None:
        val_loader = DataLoader(validation_dataset,batch_size=args.batch_size,shuffle=False,num_workers=num_workers,pin_memory=(DEVICE.type == 'cuda'))
    
    # Initialize model
    model = modelclass(args, argsGRU).to(DEVICE)
    
    # Optimizer setup
    n_losses = 5
    loss_w = torch.nn.Parameter(torch.ones(n_losses, device=DEVICE), requires_grad=True)
    optimizer = optim.AdamW(list(model.parameters()) + [loss_w], lr=args.lr, weight_decay=args.Adamw_weight_decay)
    
    # Learning rate scheduler
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=2, threshold=0.0001, threshold_mode='rel', cooldown=0, min_lr=1e-7)
    
    # TensorBoard
    writer = SummaryWriter(log_dir=f"models/logs/{args.run_name}")
    
    # Mixed precision scaler (GPU only)
    use_amp = DEVICE.type == 'cuda'
    if use_amp:
        scaler = torch.amp.GradScaler()
    
    # Gradient accumulation steps
    accum_steps = int(args.accum_steps if hasattr(args, 'accum_steps') else max(1, 8 // args.batch_size))
    print(f"Gradient accumulation steps: {accum_steps}", flush=True)
    
    # Load checkpoint if exists
    start_epoch = 0
    loss_history = []
    val_loss_history = []
    
    checkpoint_path = f"models/{args.run_name}_{MODELFILE}"
    if os.path.exists(checkpoint_path):
        checkpoint = torch.load(checkpoint_path, map_location=DEVICE, weights_only=False)
        model.load_state_dict(checkpoint['model_state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        scheduler.load_state_dict(checkpoint['scheduler_state_dict'])
        start_epoch = checkpoint['epoch']
        loss_history = checkpoint['loss']
        val_loss_history = checkpoint.get('val_loss', [])
        print(f"Model loaded from checkpoint (epoch {start_epoch}) with {count_parameters(model):,} parameters.", flush=True)
    else:
        print(f"Model initialized with {count_parameters(model):,} trainable parameters.", flush=True)
    
    # Dynamic weight adjustment setup
    if args.loss_type == 'mixed':
        #Choose 10 random indices to sample from dataset to get weights and stds
        sample_len = 10
        sample_indices = random.sample(range(0, len(dataset)), sample_len)
        sample_pairs = Subset(dataset, sample_indices)
        static_w, stds = compute_initial_weights_from_sample(sample_pairs, density_channel=0, hd_q=0.999, hist_bins=32, hist_sample_frac=0.20, device=None)
        #Each Loss is normalized to ~1
        #Boost L1 loss?
        static_w[0] = static_w[0]*10
        print(f" Initial weights and std computed with a sample of {sample_len} from dataset. ", flush=True)
        print(f"static_w: ",static_w, flush=True)
        print(f"stds: ",stds, flush=True)
        config = {
          'channel_scales': [stds[0], stds[1], stds[2], stds[3]],
          'channel_weights': [1.0, 1.0, 1.0, 1.0],
          'density_channel': 0,
          'velocity_channels': [1,2,3],
          'density_is_log': True,
          'hist_bins': 32, 'hist_vmin': -5.0, 'hist_vmax': 5.0, 'hist_sigma': None, 'hist_sample_frac': 0.2,
          'do_spectral': True, 'spec_sample_frac': 1.0,
          'hd_q': 0.999, 'hd_smooth': 1e-2
        }
        prev_losses = static_w
        prev2_losses = [1.0,1.0,1.0,1.0,1.0]
        w = static_w
    
    # Training loop
    print("=" * 80)
    print("STARTING TRAINING")
    print("=" * 80)
    training_start = time.perf_counter()
    batch_time_ema = None
    ema_alpha = 0.15
    training_phase = 1
    phase_epochs = [25, 50]
    for epoch in range(start_epoch, args.epochs):
        model.train()
        total_epoch_loss = 0.0
        sum_l1, sum_hist, sum_mass, sum_spectral, sum_hd = 0.0, 0.0, 0.0, 0.0, 0.0
        
        # Announce training phase changes
        
        if epoch in phase_epochs:
            training_phase += 1
            print("\n" + "=" * 80)
            print(f"TRAINING PHASE {training_phase} - EPOCH {epoch+1}")
            print("=" * 80 + "\n", flush=True)
            
            # Reset learning rate at phase boundaries (except first phase)
            if epoch != phase_epochs[0]:
                for g in optimizer.param_groups:
                    g['lr'] = args.lr
                print(f"Learning rate reset to {args.lr}", flush=True)
        
        # Clear cache at epoch start (GPU only)
        if DEVICE.type == 'cuda':
            torch.cuda.empty_cache()
        
        # Training batches
        optimizer.zero_grad()
        for batch_idx, batch_data in enumerate(loader):
            batch_start = time.perf_counter()
            
            # Unpack batch
            if len(batch_data) == 3:
                x, y, labels = batch_data
                if args.in_channel == 1:
                    x = x.unsqueeze(2)  # [N, t, 1, D, H, W]
                    y = y.unsqueeze(1)  # [N, 1, D, H, W]
            else:
                x, y = batch_data
            
            x, y = x.to(DEVICE, non_blocking=True), y.to(DEVICE, non_blocking=True)
            # Forward pass with mixed precision (GPU) or normal (CPU)
            if use_amp:
                with torch.amp.autocast(device_type=str(DEVICE)):
                    pred = model(x)
                # Convert to float32 for loss computation
                pred_f = pred.float()
                y_f = y.float()
            else:
                pred = model(x)
                pred_f = pred
                y_f = y
            
            # Compute losses
            l1, hist_l, spectral_l, hd_l, mass_l = return_total_loss_multichannel(y_f, pred_f, config)
            
            # Stack and weight losses
            losses = torch.stack([l1, hist_l, mass_l, spectral_l, hd_l])
            if args.loss_type == 'mixed':
                if epoch>=phase_epochs[-1]:
                    #start DWA!
                    w_t = torch.tensor(w, device=DEVICE)
                    total_loss = (w_t * losses).sum()
                else:
                    #static loss
                    total_loss = (torch.tensor(static_w, device=DEVICE) * losses).sum()
            # Scale loss for gradient accumulation
            loss = total_loss / accum_steps
            
            # Backward pass
            if use_amp:
                scaler.scale(loss).backward()
            else:
                loss.backward()
            
            # Update weights after accumulation steps
            if (batch_idx + 1) % accum_steps == 0:
                if use_amp:
                    scaler.unscale_(optimizer)
                    torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                    scaler.step(optimizer)
                    scaler.update()
                else:
                    torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                    optimizer.step()
                
                optimizer.zero_grad()
                
                # Accumulate losses
                batch_size = x.size(0)
                total_epoch_loss += total_loss.item() * batch_size
                sum_l1 += l1.item() * batch_size
                sum_hist += hist_l.item() * batch_size
                sum_mass += mass_l.item() * batch_size
                sum_spectral += spectral_l.item() * batch_size
                sum_hd += hd_l.item() * batch_size
                
                # Print progress
                batch_time_ema = print_progress(
                    batch_idx, loader, epoch, start_epoch, args,
                    total_loss, batch_start, training_start,
                    batch_time_ema, ema_alpha
                )
            
            # Periodic memory cleanup (GPU only)
            if DEVICE.type == 'cuda' and batch_idx % 5 == 0:
                torch.cuda.empty_cache()
        
        # Compute epoch averages
        N = len(dataset)
        avg_total = total_epoch_loss / N
        avg_l1 = sum_l1 / N
        avg_hist = sum_hist / N
        avg_mass = sum_mass / N
        avg_spectral = sum_spectral / N
        avg_hd = sum_hd / N
        
        # Dynamic weight adjustment
        if args.loss_type == 'mixed':
            if epoch>=phase_epochs[-1]:
                avg_losses = [avg_l1, avg_hist, avg_mass, avg_spectral, avg_hd]
                r = [prev_losses[i] / (prev2_losses[i] + 1e-8) for i in range(n_losses)]
                T = args.DWA_temperature
                K = n_losses
                exp_r = [np.exp(r_i / T) for r_i in r]
                sum_exp = sum(exp_r)
                w = [(K * e) / sum_exp for e in exp_r]
                prev2_losses = prev_losses
                prev_losses = avg_losses
                w_print = w
            else:
                w_print = static_w
        
        # Update learning rate scheduler
        last_lr = [group['lr'] for group in optimizer.param_groups]
        scheduler.step(avg_total)
        
        # Record history
        loss_history.append([avg_total, avg_l1, avg_hist, avg_mass, avg_spectral, avg_hd])
        
        # Console logging
        print(f"\nEpoch {epoch+1}/{args.epochs}: "
              f"Total={avg_total:.4f} | "
              f"L1={avg_l1:.4f} HD={avg_hd:.4f} "
              f"Mass={avg_mass:.4f} Hist={avg_hist:.4f} "
              f"Spectral={avg_spectral:.4f} | "
              f"LR={last_lr[0]:.2e} | "
              f"Weights: [{', '.join([f'{wi:.2f}' for wi in w_print])}]", flush=True)
        
        if DEVICE.type == 'cuda':
            print_memory_stats()
        
        # TensorBoard logging
        writer.add_scalar("Loss/Total", avg_total, epoch)
        writer.add_scalar("Loss/L1", avg_l1, epoch)
        writer.add_scalar("Loss/Hist", avg_hist, epoch)
        writer.add_scalar("Loss/Mass", avg_mass, epoch)
        writer.add_scalar("Loss/Spectral", avg_spectral, epoch)
        writer.add_scalar("Loss/HighDensity", avg_hd, epoch)
        writer.add_scalar("LearningRate", last_lr[0], epoch)
        
        if args.loss_type == ' mixed':
            curr_w = torch.softmax(loss_w, dim=0).detach().cpu().tolist()
            for i, name in enumerate(["L1", "Hist", "Mass", "Spectral", "HighDensity"]):
                writer.add_scalar(f"Weights/{name}", curr_w[i], epoch)
        
        # Validation every 5 epochs
        if validation_dataset is not None and epoch % 5 == 0:
            print("\nRunning validation...", flush=True)
            model.eval()
            val_total, val_l1, val_hist, val_mass, val_spectral, val_hd = 0, 0, 0, 0, 0, 0
            
            with torch.no_grad():
                for val_batch in val_loader:
                    if len(val_batch) == 3:
                        val_x, val_y, _ = val_batch
                        if args.in_channel == 1:
                            val_x = val_x.unsqueeze(2)
                            val_y = val_y.unsqueeze(1)
                    else:
                        val_x, val_y = val_batch
                    
                    val_x, val_y = val_x.to(DEVICE), val_y.to(DEVICE)
                    
                    if use_amp:
                        with torch.amp.autocast(device_type=str(DEVICE)):
                            pred = model(val_x)
                        pred = pred.float()
                        val_y = val_y.float()
                    else:
                        pred = model(val_x)
                    
                    l1, hist_l, spectral_l, hd_l, mass_l = return_total_loss_multichannel(val_y, pred, config)

                    losses = torch.stack([l1, hist_l, mass_l, spectral_l, hd_l])
                    w_t = torch.tensor(w_print, device=DEVICE)
                    total_loss = (w_t * losses).sum()
                    
                    batch_size = val_x.size(0)
                    val_total += total_loss.item() * batch_size
                    val_l1 += l1.item() * batch_size
                    val_hist += hist_l.item() * batch_size
                    val_mass += mass_l.item() * batch_size
                    val_spectral += spectral_l.item() * batch_size
                    val_hd += hd_l.item() * batch_size
            
            N_val = len(validation_dataset)
            val_loss_history.append([
                val_total / N_val, val_l1 / N_val, val_hist / N_val,
                val_mass / N_val, val_spectral / N_val, val_hd / N_val
            ])
            
            print(f"Validation: Total={val_total/N_val:.4f} | "
                  f"L1={val_l1/N_val:.4f} HD={val_hd/N_val:.4f} "
                  f"Mass={val_mass/N_val:.4f} Hist={val_hist/N_val:.4f} "
                  f"Spectral={val_spectral/N_val:.4f}\n", flush=True)
            
            # TensorBoard validation logging
            writer.add_scalar("Val_Loss/Total", val_total / N_val, epoch)
            writer.add_scalar("Val_Loss/L1", val_l1 / N_val, epoch)
            writer.add_scalar("Val_Loss/Hist", val_hist / N_val, epoch)
            writer.add_scalar("Val_Loss/Mass", val_mass / N_val, epoch)
            writer.add_scalar("Val_Loss/Spectral", val_spectral / N_val, epoch)
            writer.add_scalar("Val_Loss/HighDensity", val_hd / N_val, epoch)
        
        # Periodic checkpointing at phase boundaries
        checkpoint_epochs = [24, 49, 69, 79, 89]  # End of each phase
        if epoch in checkpoint_epochs or (epoch + 1) == args.epochs:
            print(f"\nSaving checkpoint at epoch {epoch+1}...", flush=True)
            torch.save({
                'epoch': epoch + 1,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'loss': loss_history,
                'val_loss': val_loss_history
            }, checkpoint_path)
            
            if args.loss_type == 'mixed':
                torch.save(loss_w.detach().cpu(), 
                          f"models/plots/loss_data_{args.run_name}_loss.pt")
    
    writer.close()
    
    # Final save
    print("\nSaving final model...", flush=True)
    torch.save({
        'epoch': args.epochs,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'scheduler_state_dict': scheduler.state_dict(),
        'loss': loss_history,
        'val_loss': val_loss_history
    }, checkpoint_path)
    
    # Plot training curves
    import matplotlib.pyplot as plt
    
    loss_history = np.array(loss_history)
    plt.figure(figsize=(10, 6))
    plt.plot(loss_history[:, 0], marker='o', label='Total Loss')
    plt.plot(loss_history[:, 1], marker='.', label='L1 Loss')
    plt.plot(loss_history[:, 2], marker='*', label='Histogram Loss')
    plt.plot(loss_history[:, 3], marker='+', label='Mass Loss')
    plt.plot(loss_history[:, 4], marker='x', label='Spectral Loss')
    plt.plot(loss_history[:, 5], marker='1', label='High Density Loss')
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.xlabel("Epoch")
    plt.ylabel("Avg Weighted Loss")
    plt.title(f"{args.run_name} Training Loss")
    plt.savefig(f"models/plots/training_loss_{args.run_name}.png", dpi=150)
    plt.close()
    
    # Plot validation curves if available
    if val_loss_history:
        val_loss_history = np.array(val_loss_history)
        plt.figure(figsize=(10, 6))
        plt.plot(val_loss_history[:, 0], marker='o', label='Total Loss')
        plt.plot(val_loss_history[:, 1], marker='.', label='L1 Loss')
        plt.plot(val_loss_history[:, 2], marker='*', label='Histogram Loss')
        plt.plot(val_loss_history[:, 3], marker='+', label='Mass Loss')
        plt.plot(val_loss_history[:, 4], marker='x', label='Spectral Loss')
        plt.plot(val_loss_history[:, 5], marker='1', label='High Density Loss')
        plt.yscale('log')
        plt.legend()
        plt.grid()
        plt.tight_layout()
        plt.xlabel("Validation Check (every 5 epochs)")
        plt.ylabel("Avg Weighted Loss")
        plt.title(f"{args.run_name} Validation Loss")
        plt.savefig(f"models/plots/validation_loss_{args.run_name}.png", dpi=150)
        plt.close()
    
    total_time = time.perf_counter() - training_start
    print(f"\nTraining completed in {str(datetime.timedelta(seconds=int(total_time)))}", flush=True)
    
    return min(loss_history[:, 0])  # Return minimum loss

def train_unet_unified(modelclass, dataset, args, argsGRU, validation_dataset=None):
    """
    Optimized training function with mixed precision and improved memory management and gradient accumulation.
    Notes:
    - HDF5 datasets
    - Both CPU and GPU
    - Mixed precision training (GPU only)
    - Gradient accumulation
    - Training phases with adaptive loss scaling
    - Validation
    - Learning rate scheduling
    - Periodic checkpointing
    
    Parameters:
    -----------
    dataset : torch.utils.data.Dataset
        Training dataset (HDF5Dataset)
    args : argparse.Namespace
        Training arguments
    argsGRU : argparse.Namespace
        GRU model arguments
    validation_dataset : torch.utils.data.Dataset, optional
        Validation dataset
    """
    make_dir(f"models/plots/")
    #Setup device
    DEVICE = torch.device(args.device if torch.cuda.is_available() and 'cuda' in args.device else 'cpu')
    print(f"Using device: {DEVICE}", flush=True)
    
    #Memory optimization for GPU
    if DEVICE.type == 'cuda':
        torch.cuda.empty_cache()
        torch.backends.cudnn.benchmark = False
        print_memory_stats()
    
    print(f"Training convGRUNet with {len(dataset)} samples, {args.run_name} run", flush=True)
    
    #DataLoader setup
    num_workers = args.num_workers if hasattr(args, 'num_workers') else min(4, os.cpu_count()//2)
    print(f"Using {num_workers} DataLoader workers", flush=True)
    
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, num_workers=num_workers, pin_memory=(DEVICE.type == 'cuda'), persistent_workers=(num_workers > 0))
    
    # Validation loader if provided
    if validation_dataset is not None:
        val_loader = DataLoader(validation_dataset,batch_size=args.batch_size,shuffle=False,num_workers=num_workers,pin_memory=(DEVICE.type == 'cuda'))
    
    # Initialize model
    model = modelclass(args, argsGRU).to(DEVICE)
    
    # Optimizer setup
    n_losses = 5
    loss_w = torch.nn.Parameter(torch.ones(n_losses, device=DEVICE), requires_grad=True)
    optimizer = optim.AdamW(list(model.parameters()) + [loss_w], lr=args.lr, weight_decay=args.Adamw_weight_decay)
    Earlystopper = EarlyStopping(patience=2, min_delta=0.001)
    
    # Learning rate scheduler
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=2, threshold=0.0001, threshold_mode='rel', cooldown=0, min_lr=1e-7)
    
    # TensorBoard
    writer = SummaryWriter(log_dir=f"models/logs/{args.run_name}")
    
    if hasattr(torch, 'compile'):
        print("Using torch.compile")
        model = torch.compile(
            model, 
            mode='max-autotune',  # or 'reduce-overhead' or 'default'
            backend='inductor'
        )   #First batch will be slow (compilation), then much faster
    else:
        print("torch.compile not available, upgrade to PyTorch 2.0+")

    # Mixed precision scaler (GPU only)
    use_amp = DEVICE.type == 'cuda'
    if use_amp:
        scaler = torch.amp.GradScaler()
    
    # Gradient accumulation steps
    accum_steps = int(args.accum_steps if hasattr(args, 'accum_steps') else max(1, 8 // args.batch_size))
    print(f"Gradient accumulation steps: {accum_steps}", flush=True)
    
    # Load checkpoint if exists
    start_epoch = 0
    loss_history = []
    val_loss_history = []
    
    checkpoint_path = f"models/{args.run_name}_{MODELFILE}"
    if os.path.exists(checkpoint_path):
        checkpoint = torch.load(checkpoint_path, map_location=DEVICE, weights_only=False)
        model.load_state_dict(checkpoint['model_state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        scheduler.load_state_dict(checkpoint['scheduler_state_dict'])
        start_epoch = checkpoint['epoch']
        loss_history = checkpoint['loss']
        val_loss_history = checkpoint.get('val_loss', [])
        print(f"Model loaded from checkpoint (epoch {start_epoch}) with {count_parameters(model):,} parameters.", flush=True)
    else:
        print(f"Model initialized with {count_parameters(model):,} trainable parameters.", flush=True)
    
    # Dynamic weight adjustment setup
    if args.loss_type == 'statistical':
        #Choose 50 random indices to sample from dataset to get weights and stds
        sample_len = 50
        sample_indices = random.sample(range(0, len(dataset)), sample_len)
        sample_pairs = Subset(dataset, sample_indices)
        static_w, stds = compute_initial_weights_from_sample(sample_pairs, density_channel=0, hd_q=0.999, hist_bins=32, hist_sample_frac=0.20, device=None)
        #Each Loss is normalized to ~1
        #Boost L1 loss for first 50 epochs?
        static_w[0] = static_w[0]*10
        print(f" Initial weights and std computed with a sample of {sample_len} from dataset. ", flush=True)
        print(f"static_w: ",static_w, flush=True)
        print(f"stds: ",stds, flush=True)
        config = {
          'channel_scales': [stds[0], stds[1], stds[2], stds[3]],
          'channel_weights': [1.0, 1.0, 1.0, 1.0],
          'density_channel': 0,
          'velocity_channels': [1,2,3],
          'density_is_log': True,
          'hist_bins': 32, 'hist_vmin': -5.0, 'hist_vmax': 5.0, 'hist_sigma': None, 'hist_sample_frac': 0.2,
          'do_spectral': True, 'spec_sample_frac': 1.0,
          'hd_q': 0.999, 'hd_smooth': 1e-2
        }
    
    # Training loop
    print("STARTING TRAINING")
    training_start = time.perf_counter()
    batch_time_ema = None
    ema_alpha = 0.15
    training_phase = 1
    phase_epochs = [25, 50]
    for epoch in range(start_epoch, args.epochs):
        stop_early = False
        model.train()
        total_epoch_loss = 0.0
        sum_l1, sum_hist, sum_mass, sum_spectral, sum_hd = 0.0, 0.0, 0.0, 0.0, 0.0
        
        # Announce training phase changes
        
        if epoch in phase_epochs:
            training_phase += 1
            print("\n" + "=" * 80)
            print(f"TRAINING PHASE {training_phase} - EPOCH {epoch+1}")
            print("=" * 80 + "\n", flush=True)
            
            # Reset learning rate at phase boundaries (except first phase)
            if epoch != phase_epochs[0]:
                static_w[0] = static_w[0]/10
                for g in optimizer.param_groups:
                    g['lr'] = args.lr
                print(f"Learning rate reset to {args.lr}", flush=True)
        
        # Clear cache at epoch start (GPU only)
        if DEVICE.type == 'cuda':
            torch.cuda.empty_cache()
        
        # Training batches
        optimizer.zero_grad()
        for batch_idx, batch_data in enumerate(loader):
            batch_start = time.perf_counter()
            
            # Unpack batch
            if len(batch_data) == 3:
                x, y, labels = batch_data
                if args.in_channel == 1:
                    x = x.unsqueeze(2)  # [N, t, 1, D, H, W]
                    y = y.unsqueeze(1)  # [N, 1, D, H, W]
            else:
                x, y = batch_data
            
            x, y = x.to(DEVICE, non_blocking=True), y.to(DEVICE, non_blocking=True)
            # Forward pass with mixed precision (GPU) or normal (CPU)
            if use_amp:
                with torch.amp.autocast(device_type=str(DEVICE)):
                    pred = model(x)
                # Convert to float32 for loss computation
                pred_f = pred.float()
                y_f = y.float()
            else:
                pred = model(x)
                pred_f = pred
                y_f = y
            
            # Compute losses
            l1, hist_l, spectral_l, hd_l, mass_l = return_total_loss_multichannel(y_f, pred_f, config)
            
            # Stack and weight losses
            losses = torch.stack([l1, hist_l, mass_l, spectral_l, hd_l])
            if args.loss_type == 'statistical':
                total_loss = (torch.tensor(static_w, device=DEVICE) * losses).sum()
            # Scale loss for gradient accumulation
            loss = total_loss / accum_steps
            
            # Backward pass
            if use_amp:
                scaler.scale(loss).backward()
            else:
                loss.backward()
            
            # Update weights after accumulation steps
            if (batch_idx + 1) % accum_steps == 0:
                if use_amp:
                    scaler.unscale_(optimizer)
                    torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                    scaler.step(optimizer)
                    scaler.update()
                else:
                    torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                    optimizer.step()
                
                optimizer.zero_grad()
                
                # Accumulate losses
                batch_size = x.size(0)
                total_epoch_loss += total_loss.item() * batch_size
                sum_l1 += l1.item() * batch_size
                sum_hist += hist_l.item() * batch_size
                sum_mass += mass_l.item() * batch_size
                sum_spectral += spectral_l.item() * batch_size
                sum_hd += hd_l.item() * batch_size
                
                # Print progress
                batch_time_ema = print_progress(
                    batch_idx, loader, epoch, start_epoch, args,
                    total_loss, batch_start, training_start,
                    batch_time_ema, ema_alpha
                )
            
            # Periodic memory cleanup (GPU only)
            if DEVICE.type == 'cuda' and batch_idx % 5 == 0:
                torch.cuda.empty_cache()
        
        # Compute epoch averages
        N = len(dataset)
        avg_total = total_epoch_loss / N
        avg_l1 = sum_l1 / N
        avg_hist = sum_hist / N
        avg_mass = sum_mass / N
        avg_spectral = sum_spectral / N
        avg_hd = sum_hd / N
        
        # Dynamic weight adjustment
        if args.loss_type == 'statistical':
            w_print = static_w
        
        # Update learning rate scheduler
        last_lr = [group['lr'] for group in optimizer.param_groups]
        scheduler.step(avg_total)
        
        # Record history
        loss_history.append([avg_total, avg_l1, avg_hist, avg_mass, avg_spectral, avg_hd])
        
        # Console logging
        print(f"\nEpoch {epoch+1}/{args.epochs}: "
              f"Total={avg_total:.4f} | "
              f"L1={avg_l1:.4f} HD={avg_hd:.4f} "
              f"Mass={avg_mass:.4f} Hist={avg_hist:.4f} "
              f"Spectral={avg_spectral:.4f} | "
              f"LR={last_lr[0]:.2e} | "
              f"Weights: [{', '.join([f'{wi:.2f}' for wi in w_print])}]", flush=True)
        
        if DEVICE.type == 'cuda':
            print_memory_stats()
        
        # TensorBoard logging
        writer.add_scalar("Loss/Total", avg_total, epoch)
        writer.add_scalar("Loss/L1", avg_l1, epoch)
        writer.add_scalar("Loss/Hist", avg_hist, epoch)
        writer.add_scalar("Loss/Mass", avg_mass, epoch)
        writer.add_scalar("Loss/Spectral", avg_spectral, epoch)
        writer.add_scalar("Loss/HighDensity", avg_hd, epoch)
        writer.add_scalar("LearningRate", last_lr[0], epoch)
        
        
        # Validation every 5 epochs
        if validation_dataset is not None and epoch % 5 == 0:
            print("\nRunning validation...", flush=True)
            model.eval()
            val_total, val_l1, val_hist, val_mass, val_spectral, val_hd = 0, 0, 0, 0, 0, 0
            
            with torch.no_grad():
                for val_batch in val_loader:
                    if len(val_batch) == 3:
                        val_x, val_y, _ = val_batch
                        if args.in_channel == 1:
                            val_x = val_x.unsqueeze(2)
                            val_y = val_y.unsqueeze(1)
                    else:
                        val_x, val_y = val_batch
                    
                    val_x, val_y = val_x.to(DEVICE), val_y.to(DEVICE)
                    
                    if use_amp:
                        with torch.amp.autocast(device_type=str(DEVICE)):
                            pred = model(val_x)
                        pred = pred.float()
                        val_y = val_y.float()
                    else:
                        pred = model(val_x)
                    
                    l1, hist_l, spectral_l, hd_l, mass_l = return_total_loss_multichannel(val_y, pred, config)

                    losses = torch.stack([l1, hist_l, mass_l, spectral_l, hd_l])
                    w_t = torch.tensor(w_print, device=DEVICE)
                    total_loss = (w_t * losses).sum()
                    
                    batch_size = val_x.size(0)
                    val_total += total_loss.item() * batch_size
                    val_l1 += l1.item() * batch_size
                    val_hist += hist_l.item() * batch_size
                    val_mass += mass_l.item() * batch_size
                    val_spectral += spectral_l.item() * batch_size
                    val_hd += hd_l.item() * batch_size
            
            N_val = len(validation_dataset)
            val_loss_history.append([
                val_total / N_val, val_l1 / N_val, val_hist / N_val,
                val_mass / N_val, val_spectral / N_val, val_hd / N_val
            ])
            stop_early = Earlystopper(val_total / N_val)
            
            print(f"Validation: Total={val_total/N_val:.4f} | "
                  f"L1={val_l1/N_val:.4f} HD={val_hd/N_val:.4f} "
                  f"Mass={val_mass/N_val:.4f} Hist={val_hist/N_val:.4f} "
                  f"Spectral={val_spectral/N_val:.4f}\n", flush=True)
            
            # TensorBoard validation logging
            writer.add_scalar("Val_Loss/Total", val_total / N_val, epoch)
            writer.add_scalar("Val_Loss/L1", val_l1 / N_val, epoch)
            writer.add_scalar("Val_Loss/Hist", val_hist / N_val, epoch)
            writer.add_scalar("Val_Loss/Mass", val_mass / N_val, epoch)
            writer.add_scalar("Val_Loss/Spectral", val_spectral / N_val, epoch)
            writer.add_scalar("Val_Loss/HighDensity", val_hd / N_val, epoch)
        
        # Periodic checkpointing at phase boundaries
        checkpoint_epochs = [24, 49, 69, 79, 89]  # End of each phase
        if epoch in checkpoint_epochs or (epoch + 1) == args.epochs:
            print(f"\nSaving checkpoint at epoch {epoch+1}...", flush=True)
            torch.save({
                'epoch': epoch + 1,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'loss': loss_history,
                'val_loss': val_loss_history
            }, checkpoint_path)
        
        if stop_early:
            print(f"Early stopping triggered at epoch {epoch+1}.", flush=True)
            break
                
    writer.close()
    
    # Final save
    print("\nSaving final model...", flush=True)
    torch.save({
        'epoch': args.epochs,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'scheduler_state_dict': scheduler.state_dict(),
        'loss': loss_history,
        'val_loss': val_loss_history
    }, checkpoint_path)
    
    # Plot training curves
    import matplotlib.pyplot as plt
    
    loss_history = np.array(loss_history)
    plt.figure(figsize=(10, 6))
    plt.plot(loss_history[:, 0], marker='o', label='Total Loss')
    plt.plot(loss_history[:, 1], marker='.', label='L1 Loss')
    plt.plot(loss_history[:, 2], marker='*', label='Histogram Loss')
    plt.plot(loss_history[:, 3], marker='+', label='Mass Loss')
    plt.plot(loss_history[:, 4], marker='x', label='Spectral Loss')
    plt.plot(loss_history[:, 5], marker='1', label='High Density Loss')
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.xlabel("Epoch")
    plt.ylabel("Avg Weighted Loss")
    plt.title(f"{args.run_name} Training Loss")
    plt.savefig(f"models/plots/training_loss_{args.run_name}.png", dpi=150)
    plt.close()
    
    # Plot validation curves if available
    if val_loss_history:
        val_loss_history = np.array(val_loss_history)
        plt.figure(figsize=(10, 6))
        plt.plot(val_loss_history[:, 0], marker='o', label='Total Loss')
        plt.plot(val_loss_history[:, 1], marker='.', label='L1 Loss')
        plt.plot(val_loss_history[:, 2], marker='*', label='Histogram Loss')
        plt.plot(val_loss_history[:, 3], marker='+', label='Mass Loss')
        plt.plot(val_loss_history[:, 4], marker='x', label='Spectral Loss')
        plt.plot(val_loss_history[:, 5], marker='1', label='High Density Loss')
        plt.yscale('log')
        plt.legend()
        plt.grid()
        plt.tight_layout()
        plt.xlabel("Validation Check (every 5 epochs)")
        plt.ylabel("Avg Weighted Loss")
        plt.title(f"{args.run_name} Validation Loss")
        plt.savefig(f"models/plots/validation_loss_{args.run_name}.png", dpi=150)
        plt.close()
    
    total_time = time.perf_counter() - training_start
    print(f"\nTraining completed in {str(datetime.timedelta(seconds=int(total_time)))}", flush=True)
    
    return min(loss_history[:, 0])  # Return minimum loss

def train_unet_unified_profiling(modelclass, dataset, args, argsGRU):
    """
    Optimized training function with mixed precision and improved memory management and gradient accumulation.
    Notes:
    - HDF5 datasets
    - Both CPU and GPU
    - Mixed precision training (GPU only)
    - Gradient accumulation
    - Training phases with adaptive loss scaling
    - Validation
    - Learning rate scheduling
    - Periodic checkpointing
    
    Parameters:
    -----------
    dataset : torch.utils.data.Dataset
        Training dataset (HDF5Dataset)
    args : argparse.Namespace
        Training arguments
    argsGRU : argparse.Namespace
        GRU model arguments
    validation_dataset : torch.utils.data.Dataset, optional
        Validation dataset
    """
    make_dir(f"models/plots/")
    #Setup device
    DEVICE = torch.device(args.device if torch.cuda.is_available() and 'cuda' in args.device else 'cpu')
    print(f"Using device: {DEVICE}", flush=True)
    
    #Memory optimization for GPU
    if DEVICE.type == 'cuda':
        torch.cuda.empty_cache()
        torch.backends.cudnn.benchmark = False
        print_memory_stats()
    
    print(f"Training convGRUNet with {len(dataset)} samples, {args.run_name} run", flush=True)
    
    #DataLoader setup
    num_workers = args.num_workers if hasattr(args, 'num_workers') else min(4, os.cpu_count()//2)
    print(f"Using {num_workers} DataLoader workers", flush=True)
    
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, num_workers=num_workers, pin_memory=(DEVICE.type == 'cuda'), persistent_workers=(num_workers > 0))
    
    
    # Initialize model
    model = modelclass(args, argsGRU).to(DEVICE)
    
    # Optimizer setup
    n_losses = 5
    loss_w = torch.nn.Parameter(torch.ones(n_losses, device=DEVICE), requires_grad=True)
    optimizer = optim.AdamW(list(model.parameters()) + [loss_w], lr=args.lr, weight_decay=args.Adamw_weight_decay)
    
    # Learning rate scheduler
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=2, threshold=0.0001, threshold_mode='rel', cooldown=0, min_lr=1e-7)
    
    # TensorBoard
    writer = SummaryWriter(log_dir=f"models/logs/{args.run_name}")
    
    # Mixed precision scaler (GPU only)
    use_amp = DEVICE.type == 'cuda'
    if use_amp:
        scaler = torch.amp.GradScaler()
    
    # Gradient accumulation steps
    accum_steps = int(args.accum_steps if hasattr(args, 'accum_steps') else max(1, 8 // args.batch_size))
    print(f"Gradient accumulation steps: {accum_steps}", flush=True)
    
    # Load checkpoint if exists
    start_epoch = 0
    loss_history = []

    print(f"Model initialized with {count_parameters(model):,} trainable parameters.", flush=True)
    
    # Dynamic weight adjustment setup
    if args.loss_type == 'statistical':
        #Choose 50 random indices to sample from dataset to get weights and stds
        sample_len = 50
        sample_indices = random.sample(range(0, len(dataset)), sample_len)
        sample_pairs = Subset(dataset, sample_indices)
        static_w, stds = compute_initial_weights_from_sample(sample_pairs, density_channel=0, hd_q=0.999, hist_bins=32, hist_sample_frac=0.20, device=None)
        #Each Loss is normalized to ~1
        #Boost L1 loss for first 50 epochs?
        static_w[0] = static_w[0]*10
        print(f" Initial weights and std computed with a sample of {sample_len} from dataset. ", flush=True)
        print(f"static_w: ",static_w, flush=True)
        print(f"stds: ",stds, flush=True)
        config = {
          'channel_scales': [stds[0], stds[1], stds[2], stds[3]],
          'channel_weights': [1.0, 1.0, 1.0, 1.0],
          'density_channel': 0,
          'velocity_channels': [1,2,3],
          'density_is_log': True,
          'hist_bins': 32, 'hist_vmin': -5.0, 'hist_vmax': 5.0, 'hist_sigma': None, 'hist_sample_frac': 0.2,
          'do_spectral': True, 'spec_sample_frac': 1.0,
          'hd_q': 0.999, 'hd_smooth': 1e-2
        }
    
    # Training loop
    print("=" * 80)
    print("STARTING TRAINING")
    print("=" * 80)
    training_start = time.perf_counter()
    batch_time_ema = None
    ema_alpha = 0.15
    for epoch in range(start_epoch, args.epochs):
        model.train()
        total_epoch_loss = 0.0
        sum_l1, sum_hist, sum_mass, sum_spectral, sum_hd = 0.0, 0.0, 0.0, 0.0, 0.0
        
        # Clear cache at epoch start (GPU only)
        if DEVICE.type == 'cuda':
            torch.cuda.empty_cache()
        
        # Training batches
        optimizer.zero_grad()
        for batch_idx, batch_data in enumerate(loader):
            batch_start = time.perf_counter()
            
            # Unpack batch
            if len(batch_data) == 3:
                x, y, labels = batch_data
                if args.in_channel == 1:
                    x = x.unsqueeze(2)  # [N, t, 1, D, H, W]
                    y = y.unsqueeze(1)  # [N, 1, D, H, W]
            else:
                x, y = batch_data
            
            x, y = x.to(DEVICE, non_blocking=True), y.to(DEVICE, non_blocking=True)
            # Forward pass with mixed precision (GPU) or normal (CPU)
            if use_amp:
                with torch.amp.autocast(device_type=str(DEVICE)):
                    pred = model(x)
                # Convert to float32 for loss computation
                pred_f = pred.float()
                y_f = y.float()
            else:
                pred = model(x)
                pred_f = pred
                y_f = y
            
            # Compute losses
            l1, hist_l, spectral_l, hd_l, mass_l = return_total_loss_multichannel(y_f, pred_f, config)
            
            # Stack and weight losses
            losses = torch.stack([l1, hist_l, mass_l, spectral_l, hd_l])
            if args.loss_type == 'statistical':
                total_loss = (torch.tensor(static_w, device=DEVICE) * losses).sum()
            # Scale loss for gradient accumulation
            loss = total_loss / accum_steps
            
            # Backward pass
            if use_amp:
                scaler.scale(loss).backward()
            else:
                loss.backward()
            
            # Update weights after accumulation steps
            if (batch_idx + 1) % accum_steps == 0:
                if use_amp:
                    scaler.unscale_(optimizer)
                    torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                    scaler.step(optimizer)
                    scaler.update()
                else:
                    torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                    optimizer.step()
                
                optimizer.zero_grad()
                
                # Accumulate losses
                batch_size = x.size(0)
                total_epoch_loss += total_loss.item() * batch_size
                sum_l1 += l1.item() * batch_size
                sum_hist += hist_l.item() * batch_size
                sum_mass += mass_l.item() * batch_size
                sum_spectral += spectral_l.item() * batch_size
                sum_hd += hd_l.item() * batch_size
                
                # Print progress
                batch_time_ema = print_progress(
                    batch_idx, loader, epoch, start_epoch, args,
                    total_loss, batch_start, training_start,
                    batch_time_ema, ema_alpha
                )
            
            # Periodic memory cleanup (GPU only)
            if DEVICE.type == 'cuda' and batch_idx % 5 == 0:
                torch.cuda.empty_cache()
        
        # Compute epoch averages
        N = len(dataset)
        avg_total = total_epoch_loss / N
        avg_l1 = sum_l1 / N
        avg_hist = sum_hist / N
        avg_mass = sum_mass / N
        avg_spectral = sum_spectral / N
        avg_hd = sum_hd / N
        
        # Dynamic weight adjustment
        if args.loss_type == 'statistical':
            w_print = static_w
        
        # Update learning rate scheduler
        last_lr = [group['lr'] for group in optimizer.param_groups]
        scheduler.step(avg_total)
        
        # Record history
        loss_history.append([avg_total, avg_l1, avg_hist, avg_mass, avg_spectral, avg_hd])
        
        # Console logging
        print(f"\nEpoch {epoch+1}/{args.epochs}: "
              f"Total={avg_total:.4f} | "
              f"L1={avg_l1:.4f} HD={avg_hd:.4f} "
              f"Mass={avg_mass:.4f} Hist={avg_hist:.4f} "
              f"Spectral={avg_spectral:.4f} | "
              f"LR={last_lr[0]:.2e} | "
              f"Weights: [{', '.join([f'{wi:.2f}' for wi in w_print])}]", flush=True)
        
        if DEVICE.type == 'cuda':
            print_memory_stats()
        
        # TensorBoard logging
        writer.add_scalar("Loss/Total", avg_total, epoch)
        writer.add_scalar("Loss/L1", avg_l1, epoch)
        writer.add_scalar("Loss/Hist", avg_hist, epoch)
        writer.add_scalar("Loss/Mass", avg_mass, epoch)
        writer.add_scalar("Loss/Spectral", avg_spectral, epoch)
        writer.add_scalar("Loss/HighDensity", avg_hd, epoch)
        writer.add_scalar("LearningRate", last_lr[0], epoch)
        
    writer.close()

    total_time = time.perf_counter() - training_start
    print(f"\nTraining completed in {str(datetime.timedelta(seconds=int(total_time)))}", flush=True)
    
    return loss_history

def train_unet_unified_delta(modelclass, dataset, args, argsGRU, validation_dataset=None, dont_save=False):
    """
    Optimized training function with delta training with mixed precision and improved memory management and gradient accumulation.
    Notes:
    - HDF5 datasets
    - Both CPU and GPU
    - Mixed precision training (GPU only)
    - Gradient accumulation
    - Training phases with adaptive loss scaling
    - Validation
    - Learning rate scheduling
    - Periodic checkpointing
    
    Parameters:
    -----------
    dataset : torch.utils.data.Dataset
        Training dataset (HDF5Dataset)
    args : argparse.Namespace
        Training arguments
    argsGRU : argparse.Namespace
        GRU model arguments
    validation_dataset : torch.utils.data.Dataset, optional
        Validation dataset
    """
    make_dir(f"models/plots/")
    #Setup device
    DEVICE = torch.device(args.device if torch.cuda.is_available() and 'cuda' in args.device else 'cpu')
    print(f"Using device: {DEVICE}", flush=True)
    
    #Memory optimization for GPU
    if DEVICE.type == 'cuda':
        torch.cuda.empty_cache()
        torch.backends.cudnn.benchmark = False
        print_memory_stats()
    
    print(f"Training convGRUNet with {len(dataset)} samples, {args.run_name} run", flush=True)
    
    #DataLoader setup
    num_workers = args.num_workers if hasattr(args, 'num_workers') else min(4, os.cpu_count()//2)
    print(f"Using {num_workers} DataLoader workers", flush=True)
    
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, num_workers=num_workers, pin_memory=(DEVICE.type == 'cuda'), persistent_workers=(num_workers > 0))
    
    # Validation loader if provided
    if validation_dataset is not None:
        val_loader = DataLoader(validation_dataset,batch_size=args.batch_size,shuffle=False,num_workers=num_workers,pin_memory=(DEVICE.type == 'cuda'))
    
    # Initialize model
    model = modelclass(args, argsGRU).to(DEVICE)
    
    # Optimizer setup
    n_losses = 5
    loss_w = torch.nn.Parameter(torch.ones(n_losses, device=DEVICE), requires_grad=True)
    optimizer = optim.AdamW(list(model.parameters()) + [loss_w], lr=args.lr, weight_decay=args.Adamw_weight_decay)
    Earlystopper = EarlyStopping(patience=2, min_delta=0.001)
    
    # Learning rate scheduler
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=2, threshold=0.0001, threshold_mode='rel', cooldown=0, min_lr=1e-7)
    
    # TensorBoard
    writer = SummaryWriter(log_dir=f"models/logs/{args.run_name}")
    
    if hasattr(torch, 'compile'):
        print("Using torch.compile")
        model = torch.compile(
            model, 
            mode='max-autotune',  # or 'reduce-overhead' or 'default'
            backend='inductor'
        )   #First batch will be slow (compilation), then much faster
    else:
        print("torch.compile not available, upgrade to PyTorch 2.0+")
            
    # Mixed precision scaler (GPU only)
    use_amp = DEVICE.type == 'cuda'
    if use_amp:
        scaler = torch.amp.GradScaler()
    
    # Gradient accumulation steps
    accum_steps = int(args.accum_steps if hasattr(args, 'accum_steps') else max(1, 8 // args.batch_size))
    print(f"Gradient accumulation steps: {accum_steps}", flush=True)
    
    # Load checkpoint if exists
    start_epoch = 0
    loss_history = []
    val_loss_history = []
    
    checkpoint_path = f"models/{args.run_name}_{MODELFILE}"
    if os.path.exists(checkpoint_path):
        checkpoint = torch.load(checkpoint_path, map_location=DEVICE, weights_only=False)
        model.load_state_dict(checkpoint['model_state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        scheduler.load_state_dict(checkpoint['scheduler_state_dict'])
        start_epoch = checkpoint['epoch']
        loss_history = checkpoint['loss']
        val_loss_history = checkpoint.get('val_loss', [])
        print(f"Model loaded from checkpoint (epoch {start_epoch}) with {count_parameters(model):,} parameters.", flush=True)
    else:
        print(f"Model initialized with {count_parameters(model):,} trainable parameters.", flush=True)
    
    # Dynamic weight adjustment setup
    if args.loss_type == 'statistical':
        #Choose 50 random indices to sample from dataset to get weights and stds
        sample_len = 50
        sample_indices = random.sample(range(0, len(dataset)), sample_len)
        sample_pairs = Subset(dataset, sample_indices)
        static_w, stds = compute_initial_weights_from_sample(sample_pairs, density_channel=0, hd_q=0.999, hist_bins=32, hist_sample_frac=0.20, device=None)
        #Each Loss is normalized to ~1
        #Boost L1 loss for first 50 epochs?
        static_w[0] = static_w[0]*10
        print(f" Initial weights and std computed with a sample of {sample_len} from dataset. ", flush=True)
        print(f"static_w: ",static_w, flush=True)
        print(f"stds: ",stds, flush=True)
        config = {
          'type':'delta',
          'channel_scales': [stds[0], stds[1], stds[2], stds[3]],
          'channel_weights': [1.0, 1.0, 1.0, 1.0],
          'density_channel': 0,
          'velocity_channels': [1,2,3],
          'density_is_log': True,
          'hist_bins': 32, 'hist_vmin': -5.0, 'hist_vmax': 5.0, 'hist_sigma': None, 'hist_sample_frac': 0.2,
          'do_spectral': True, 'spec_sample_frac': 1.0,
          'hd_q': 0.999, 'hd_smooth': 1e-2
        }
    
    # Training loop
    print("STARTING TRAINING")
    training_start = time.perf_counter()
    batch_time_ema = None
    ema_alpha = 0.15
    training_phase = 1
    phase_epochs = [25, 50]
    for epoch in range(start_epoch, args.epochs):
        stop_early = False
        model.train()
        total_epoch_loss = 0.0
        sum_l1, sum_hist, sum_mass, sum_spectral, sum_hd = 0.0, 0.0, 0.0, 0.0, 0.0
        
        # Announce training phase changes
        
        if epoch in phase_epochs:
            training_phase += 1
            print("\n" + "=" * 80)
            print(f"TRAINING PHASE {training_phase} - EPOCH {epoch+1}")
            print("=" * 80 + "\n", flush=True)
            
            # Reset learning rate at phase boundaries (except first phase)
            if epoch != phase_epochs[0]:
                static_w[0] = static_w[0]/10
                for g in optimizer.param_groups:
                    g['lr'] = args.lr
                print(f"Learning rate reset to {args.lr}", flush=True)
        
        # Clear cache at epoch start (GPU only)
        if DEVICE.type == 'cuda':
            torch.cuda.empty_cache()
        
        # Training batches
        optimizer.zero_grad()
        for batch_idx, batch_data in enumerate(loader):
            batch_start = time.perf_counter()
            
            # Unpack batch
            if len(batch_data) == 3:
                x, y, labels = batch_data
                if args.in_channel == 1:
                    x = x.unsqueeze(2)  # [N, t, 1, D, H, W]
                    y = y.unsqueeze(1)  # [N, 1, D, H, W]
            else:
                x, y = batch_data
            
            x, y = x.to(DEVICE, non_blocking=True), y.to(DEVICE, non_blocking=True)
            delta = y-x
            # Forward pass with mixed precision (GPU) or normal (CPU)
            if use_amp:
                with torch.amp.autocast(device_type=str(DEVICE)):
                    deltapred = model(x)
                # Convert to float32 for loss computation
                pred_f = deltapred.float()
                delta_f = delta.float()
            else:
                deltapred = model(x)
                pred_f = deltapred
                delta_f = delta
            
            # Compute losses
            l1, hist_l, spectral_l, hd_l, mass_l = return_total_loss_multichannel(delta_f, pred_f, config, input_state=x)
            
            # Stack and weight losses
            losses = torch.stack([l1, hist_l, mass_l, spectral_l, hd_l])
            if args.loss_type == 'statistical':
                total_loss = (torch.tensor(static_w, device=DEVICE) * losses).sum()
            # Scale loss for gradient accumulation
            loss = total_loss / accum_steps
            
            # Backward pass
            if use_amp:
                scaler.scale(loss).backward()
            else:
                loss.backward()
            
            # Update weights after accumulation steps
            if (batch_idx + 1) % accum_steps == 0:
                if use_amp:
                    scaler.unscale_(optimizer)
                    torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                    scaler.step(optimizer)
                    scaler.update()
                else:
                    torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                    optimizer.step()
                
                optimizer.zero_grad()
                
                # Accumulate losses
                batch_size = x.size(0)
                total_epoch_loss += total_loss.item() * batch_size
                sum_l1 += l1.item() * batch_size
                sum_hist += hist_l.item() * batch_size
                sum_mass += mass_l.item() * batch_size
                sum_spectral += spectral_l.item() * batch_size
                sum_hd += hd_l.item() * batch_size
                
                # Print progress
                batch_time_ema = print_progress(
                    batch_idx, loader, epoch, start_epoch, args,
                    total_loss, batch_start, training_start,
                    batch_time_ema, ema_alpha
                )
            
            # Periodic memory cleanup (GPU only)
            if DEVICE.type == 'cuda' and batch_idx % 5 == 0:
                torch.cuda.empty_cache()
        
        # Compute epoch averages
        N = len(dataset)
        avg_total = total_epoch_loss / N
        avg_l1 = sum_l1 / N
        avg_hist = sum_hist / N
        avg_mass = sum_mass / N
        avg_spectral = sum_spectral / N
        avg_hd = sum_hd / N
        
        # Dynamic weight adjustment
        if args.loss_type == 'statistical':
            w_print = static_w
        
        # Update learning rate scheduler
        last_lr = [group['lr'] for group in optimizer.param_groups]
        scheduler.step(avg_total)
        
        # Record history
        loss_history.append([avg_total, avg_l1, avg_hist, avg_mass, avg_spectral, avg_hd])
        
        # Console logging
        print(f"\nEpoch {epoch+1}/{args.epochs}: "
              f"Total={avg_total:.4f} | "
              f"L1={avg_l1:.4f} HD={avg_hd:.4f} "
              f"Mass={avg_mass:.4f} Hist={avg_hist:.4f} "
              f"Spectral={avg_spectral:.4f} | "
              f"LR={last_lr[0]:.2e} | "
              f"Weights: [{', '.join([f'{wi:.2f}' for wi in w_print])}]", flush=True)
        
        if DEVICE.type == 'cuda':
            print_memory_stats()
        
        # TensorBoard logging
        writer.add_scalar("Loss/Total", avg_total, epoch)
        writer.add_scalar("Loss/L1", avg_l1, epoch)
        writer.add_scalar("Loss/Hist", avg_hist, epoch)
        writer.add_scalar("Loss/Mass", avg_mass, epoch)
        writer.add_scalar("Loss/Spectral", avg_spectral, epoch)
        writer.add_scalar("Loss/HighDensity", avg_hd, epoch)
        writer.add_scalar("LearningRate", last_lr[0], epoch)
        
        
        # Validation every 5 epochs
        if validation_dataset is not None and epoch % 5 == 0:
            print("\nRunning validation...", flush=True)
            model.eval()
            val_total, val_l1, val_hist, val_mass, val_spectral, val_hd = 0, 0, 0, 0, 0, 0
            
            with torch.no_grad():
                for val_batch in val_loader:
                    if len(val_batch) == 3:
                        val_x, val_y, _ = val_batch
                        if args.in_channel == 1:
                            val_x = val_x.unsqueeze(2)
                            val_y = val_y.unsqueeze(1)
                    else:
                        val_x, val_y = val_batch
                    
                    val_x, val_y = val_x.to(DEVICE), val_y.to(DEVICE)
                    val_delta = val_y-val_x
                    if use_amp:
                        with torch.amp.autocast(device_type=str(DEVICE)):
                            val_pred = model(val_x)
                        pred = val_pred.float()
                        val_delta_f = val_delta.float()
                    else:
                        pred = model(val_x)
                        val_delta_f = val_delta
                    
                    l1, hist_l, spectral_l, hd_l, mass_l = return_total_loss_multichannel(val_delta_f, pred, config, input_state=val_x)

                    losses = torch.stack([l1, hist_l, mass_l, spectral_l, hd_l])
                    w_t = torch.tensor(w_print, device=DEVICE)
                    total_loss = (w_t * losses).sum()
                    
                    batch_size = val_x.size(0)
                    val_total += total_loss.item() * batch_size
                    val_l1 += l1.item() * batch_size
                    val_hist += hist_l.item() * batch_size
                    val_mass += mass_l.item() * batch_size
                    val_spectral += spectral_l.item() * batch_size
                    val_hd += hd_l.item() * batch_size
            
            N_val = len(validation_dataset)
            val_loss_history.append([
                val_total / N_val, val_l1 / N_val, val_hist / N_val,
                val_mass / N_val, val_spectral / N_val, val_hd / N_val
            ])
            stop_early = Earlystopper(val_total / N_val)
            
            print(f"Validation: Total={val_total/N_val:.4f} | "
                  f"L1={val_l1/N_val:.4f} HD={val_hd/N_val:.4f} "
                  f"Mass={val_mass/N_val:.4f} Hist={val_hist/N_val:.4f} "
                  f"Spectral={val_spectral/N_val:.4f}\n", flush=True)
            
            # TensorBoard validation logging
            writer.add_scalar("Val_Loss/Total", val_total / N_val, epoch)
            writer.add_scalar("Val_Loss/L1", val_l1 / N_val, epoch)
            writer.add_scalar("Val_Loss/Hist", val_hist / N_val, epoch)
            writer.add_scalar("Val_Loss/Mass", val_mass / N_val, epoch)
            writer.add_scalar("Val_Loss/Spectral", val_spectral / N_val, epoch)
            writer.add_scalar("Val_Loss/HighDensity", val_hd / N_val, epoch)
        
        # Periodic checkpointing at phase boundaries
        checkpoint_epochs = [24, 49, 69, 79, 89]  # End of each phase
        if epoch in checkpoint_epochs:
            print(f"\nSaving checkpoint at epoch {epoch+1}...", flush=True)
            torch.save({
                'epoch': epoch + 1,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'loss': loss_history,
                'val_loss': val_loss_history
            }, checkpoint_path)
        
        if stop_early:
            print(f"Early stopping triggered at epoch {epoch+1}.", flush=True)
            break
                
    writer.close()
    loss_history = np.array(loss_history)
    if dont_save:
        return min(loss_history[:, 0])  # Return minimum loss
    # Final save
    print("\nSaving final model...", flush=True)
    torch.save({
        'epoch': args.epochs,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'scheduler_state_dict': scheduler.state_dict(),
        'loss': loss_history,
        'val_loss': val_loss_history
    }, checkpoint_path)
    
    # Plot training curves
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.plot(loss_history[:, 0], marker='o', label='Total Loss')
    plt.plot(loss_history[:, 1], marker='.', label='L1 Loss')
    plt.plot(loss_history[:, 2], marker='*', label='Histogram Loss')
    plt.plot(loss_history[:, 3], marker='+', label='Mass Loss')
    plt.plot(loss_history[:, 4], marker='x', label='Spectral Loss')
    plt.plot(loss_history[:, 5], marker='1', label='High Density Loss')
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.xlabel("Epoch")
    plt.ylabel("Avg Weighted Loss")
    plt.title(f"{args.run_name} Training Loss")
    plt.savefig(f"models/plots/training_loss_{args.run_name}.png", dpi=150)
    plt.close()
    
    # Plot validation curves if available
    if val_loss_history:
        val_loss_history = np.array(val_loss_history)
        plt.figure(figsize=(10, 6))
        plt.plot(val_loss_history[:, 0], marker='o', label='Total Loss')
        plt.plot(val_loss_history[:, 1], marker='.', label='L1 Loss')
        plt.plot(val_loss_history[:, 2], marker='*', label='Histogram Loss')
        plt.plot(val_loss_history[:, 3], marker='+', label='Mass Loss')
        plt.plot(val_loss_history[:, 4], marker='x', label='Spectral Loss')
        plt.plot(val_loss_history[:, 5], marker='1', label='High Density Loss')
        plt.yscale('log')
        plt.legend()
        plt.grid()
        plt.tight_layout()
        plt.xlabel("Validation Check (every 5 epochs)")
        plt.ylabel("Avg Weighted Loss")
        plt.title(f"{args.run_name} Validation Loss")
        plt.savefig(f"models/plots/validation_loss_{args.run_name}.png", dpi=150)
        plt.close()
    
    total_time = time.perf_counter() - training_start
    print(f"\nTraining completed in {str(datetime.timedelta(seconds=int(total_time)))}", flush=True)
    
    return min(loss_history[:, 0])  # Return minimum loss

def train_unet_unified_delta_chunked(modelclass, chunk_manager, args, argsGRU, validate_dataset=True, dont_save=False):
    """
    Optimized training function with delta training with mixed precision and improved memory management and gradient accumulation.
    Notes:
    - HDF5 datasets
    - Both CPU and GPU
    - Mixed precision training (GPU only)
    - Gradient accumulation
    - Training phases with adaptive loss scaling
    - Validation
    - Learning rate scheduling
    - Periodic checkpointing
    
    Parameters:
    -----------
    dataset : torch.utils.data.Dataset
        Training dataset (HDF5Dataset)
    args : argparse.Namespace
        Training arguments
    argsGRU : argparse.Namespace
        GRU model arguments
    validation_dataset : torch.utils.data.Dataset, optional
        Validation dataset
    """
    make_dir(f"models/plots/")
    #Setup device
    DEVICE = torch.device(args.device if torch.cuda.is_available() and 'cuda' in args.device else 'cpu')
    print(f"Using device: {DEVICE}", flush=True)
    
    #Memory optimization for GPU
    if DEVICE.type == 'cuda':
        torch.cuda.empty_cache()
        torch.backends.cudnn.benchmark = False
        print_memory_stats()
    
    print(f"Training convGRUNet with {chunk_manager.total_samples} samples, {args.run_name} run", flush=True)
    
    # Validation loader if provided
    if validate_dataset:
        print("Loading Validation Set into RAM...", flush=True)
        val_dataset = chunk_manager.load_data_to_ram(chunk_manager.test_keys, desc="Validation")
        val_loader = DataLoader(val_dataset, batch_size=args.batch_size, shuffle=False, num_workers=0,pin_memory=(DEVICE.type == 'cuda'))
    
    # Initialize model
    model = modelclass(args, argsGRU).to(DEVICE)
    
    # Optimizer setup
    n_losses = 5
    loss_w = torch.nn.Parameter(torch.ones(n_losses, device=DEVICE), requires_grad=True)
    optimizer = optim.AdamW(list(model.parameters()) + [loss_w], lr=args.lr, weight_decay=args.Adamw_weight_decay)
    Earlystopper = EarlyStopping(patience=2, min_delta=0.001)
    
    # Learning rate scheduler
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=2, threshold=0.0001, threshold_mode='rel', cooldown=0, min_lr=1e-7)
    
    # TensorBoard
    writer = SummaryWriter(log_dir=f"models/logs/{args.run_name}")
    
    if hasattr(torch, 'compile'):
        print("Using torch.compile")
        model = torch.compile(
            model, 
            mode='max-autotune',  # or 'reduce-overhead' or 'default'
            backend='inductor'
        )   #First batch will be slow (compilation), then much faster
    else:
        print("torch.compile not available, upgrade to PyTorch 2.0+")
            
    # Mixed precision scaler (GPU only)
    use_amp = DEVICE.type == 'cuda'
    if use_amp:
        scaler = torch.amp.GradScaler()
    
    # Gradient accumulation steps
    accum_steps = int(args.accum_steps if hasattr(args, 'accum_steps') else max(1, 8 // args.batch_size))
    print(f"Gradient accumulation steps: {accum_steps}", flush=True)
    
    # Load checkpoint if exists
    start_epoch = 0
    loss_history = []
    val_loss_history = []

    checkpoint_path = f"models/{args.run_name}_{MODELFILE}"
    if os.path.exists(checkpoint_path):
        checkpoint = torch.load(checkpoint_path, map_location=DEVICE, weights_only=False)
        model.load_state_dict(checkpoint['model_state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        scheduler.load_state_dict(checkpoint['scheduler_state_dict'])
        start_epoch = checkpoint['epoch']
        loss_history = checkpoint['loss']
        val_loss_history = checkpoint.get('val_loss', [])
        print(f"Model loaded from checkpoint (epoch {start_epoch}) with {count_parameters(model):,} parameters.", flush=True)
    else:
        print(f"Model initialized with {count_parameters(model):,} trainable parameters.", flush=True)
    
    # Dynamic weight adjustment setup
    first_chunk_dataset = None
    if args.loss_type == 'statistical':
        #Calculate Statistics (Using ONLY the first training chunk to save time)
        print("Loading First Chunk for Statistics...", flush=True)
        first_chunk_dataset = chunk_manager.load_data_to_ram(chunk_manager.train_chunks[0], desc="Chunk 0")
        print("Computing initial weights from first chunk...", flush=True)
        sample_pairs = []
        #Choose 50 random indices to sample from dataset to get weights and stds
        sample_len = 50
        for i in range(min(len(first_chunk_dataset), sample_len)): # Use subset for speed
            sample_pairs.append(first_chunk_dataset[i])
            
        static_w, stds = compute_initial_weights_from_sample(sample_pairs, density_channel=0, hd_q=0.999, hist_bins=32, hist_sample_frac=0.20, device=DEVICE)
        #Each Loss is normalized to ~1
        #Boost L1 loss for first 50 epochs?
        static_w[0] = static_w[0]*10
        print(f" Initial weights and std computed with a sample of {sample_len} from dataset. ", flush=True)
        print(f"static_w: ",static_w, flush=True)
        print(f"stds: ",stds, flush=True)
        config = {
          'type':'delta',
          'channel_scales': [stds[0], stds[1], stds[2], stds[3]],
          'channel_weights': [1.0, 1.0, 1.0, 1.0],
          'density_channel': 0,
          'velocity_channels': [1,2,3],
          'density_is_log': True,
          'hist_bins': 32, 'hist_vmin': -5.0, 'hist_vmax': 5.0, 'hist_sigma': None, 'hist_sample_frac': 0.2,
          'do_spectral': True, 'spec_sample_frac': 1.0,
          'hd_q': 0.999, 'hd_smooth': 1e-2
        }
        del sample_pairs
    
    # Training loop
    print(f"Starting Training on {len(chunk_manager.train_chunks)} chunks per epoch...", flush=True)
    training_start = time.perf_counter()
    batch_time_ema = None
    ema_alpha = 0.15
    training_phase = 1
    phase_epochs = [25, 50]
    for epoch in range(start_epoch, args.epochs):
        stop_early = False
        model.train()
        total_epoch_loss = 0.0
        sum_l1, sum_hist, sum_mass, sum_spectral, sum_hd = 0.0, 0.0, 0.0, 0.0, 0.0
        
        # Announce training phase changes
        if epoch in phase_epochs:
            training_phase += 1
            print("\n" + "=" * 80)
            print(f"TRAINING PHASE {training_phase} - EPOCH {epoch+1}")
            print("=" * 80 + "\n", flush=True)
            
            # Reset learning rate at phase boundaries (except first phase)
            if epoch != phase_epochs[0]:
                static_w[0] = static_w[0]/10
                for g in optimizer.param_groups:
                    g['lr'] = args.lr
                print(f"Learning rate reset to {args.lr}", flush=True)
        
        # Clear cache at epoch start (GPU only)
        if DEVICE.type == 'cuda':
            torch.cuda.empty_cache()
        N=0
        # Iterate over CHUNKS
        for i, chunk_keys in enumerate(chunk_manager.train_chunks):
            if first_chunk_dataset is not None and i==0:
                train_dataset = first_chunk_dataset
            else:
                train_dataset = chunk_manager.load_data_to_ram(chunk_keys, desc=f"Ep {epoch} - Chunk {i}")
            train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True, num_workers=0, pin_memory=(DEVICE.type == 'cuda'))
            # Training batches
            for batch_idx, batch_data in enumerate(train_loader):
                batch_start = time.perf_counter()
                
                # Unpack batch
                if len(batch_data) == 3:
                    x, y, labels = batch_data
                    if args.in_channel == 1:
                        x = x.unsqueeze(2)  # [N, t, 1, D, H, W]
                        y = y.unsqueeze(1)  # [N, 1, D, H, W]
                else:
                    x, y = batch_data
                
                x, y = x.to(DEVICE, non_blocking=True), y.to(DEVICE, non_blocking=True)
                delta = y-x
                # Forward pass with mixed precision (GPU) or normal (CPU)
                optimizer.zero_grad()
                if use_amp:
                    with torch.amp.autocast(device_type=str(DEVICE)):
                        deltapred = model(x)
                    # Convert to float32 for loss computation
                    pred_f = deltapred.float()
                    delta_f = delta.float()
                else:
                    deltapred = model(x)
                    pred_f = deltapred
                    delta_f = delta
                
                # Compute losses
                l1, hist_l, spectral_l, hd_l, mass_l = return_total_loss_multichannel(delta_f, pred_f, config, input_state=x)
                
                # Stack and weight losses
                losses = torch.stack([l1, hist_l, mass_l, spectral_l, hd_l])
                if args.loss_type == 'statistical':
                    total_loss = (torch.tensor(static_w, device=DEVICE) * losses).sum()
                # Scale loss for gradient accumulation
                loss = total_loss / accum_steps
                
                # Backward pass
                if use_amp:
                    scaler.scale(loss).backward()
                else:
                    loss.backward()
                
                # Update weights after accumulation steps
                if (batch_idx + 1) % accum_steps == 0:
                    if use_amp:
                        scaler.unscale_(optimizer)
                        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                        scaler.step(optimizer)
                        scaler.update()
                    else:
                        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                        optimizer.step()
                    
                    optimizer.zero_grad()
                    
                    # Accumulate losses
                    batch_size = x.size(0)
                    total_epoch_loss += total_loss.item() * batch_size
                    sum_l1 += l1.item() * batch_size
                    sum_hist += hist_l.item() * batch_size
                    sum_mass += mass_l.item() * batch_size
                    sum_spectral += spectral_l.item() * batch_size
                    sum_hd += hd_l.item() * batch_size
                    
                    # Print progress
                    batch_time_ema = print_progress(
                        batch_idx, train_loader, epoch, start_epoch, args,
                        total_loss, batch_start, training_start,
                        batch_time_ema, ema_alpha
                    )
                
                # Periodic memory cleanup (GPU only)
                if DEVICE.type == 'cuda' and batch_idx % 5 == 0:
                    torch.cuda.empty_cache()
            
            N += len(train_dataset)
            #Delete dataset to free RAM for next chunk
            del train_dataset
            del train_loader
            
        # Compute epoch averages
        avg_total = total_epoch_loss / N
        avg_l1 = sum_l1 / N
        avg_hist = sum_hist / N
        avg_mass = sum_mass / N
        avg_spectral = sum_spectral / N
        avg_hd = sum_hd / N
        
        # Dynamic weight adjustment
        if args.loss_type == 'statistical':
            w_print = static_w
        
        # Update learning rate scheduler
        last_lr = [group['lr'] for group in optimizer.param_groups]
        scheduler.step(avg_total)
        
        # Record history
        loss_history.append([avg_total, avg_l1, avg_hist, avg_mass, avg_spectral, avg_hd])
        
        # Console logging
        print(f"\nEpoch {epoch+1}/{args.epochs}: "
              f"Total={avg_total:.4f} | "
              f"L1={avg_l1:.4f} HD={avg_hd:.4f} "
              f"Mass={avg_mass:.4f} Hist={avg_hist:.4f} "
              f"Spectral={avg_spectral:.4f} | "
              f"LR={last_lr[0]:.2e} | "
              f"Weights: [{', '.join([f'{wi:.2f}' for wi in w_print])}]", flush=True)
        
        if DEVICE.type == 'cuda':
            print_memory_stats()
        
        # TensorBoard logging
        writer.add_scalar("Loss/Total", avg_total, epoch)
        writer.add_scalar("Loss/L1", avg_l1, epoch)
        writer.add_scalar("Loss/Hist", avg_hist, epoch)
        writer.add_scalar("Loss/Mass", avg_mass, epoch)
        writer.add_scalar("Loss/Spectral", avg_spectral, epoch)
        writer.add_scalar("Loss/HighDensity", avg_hd, epoch)
        writer.add_scalar("LearningRate", last_lr[0], epoch)
        
        
        # Validation every 5 epochs
        if validate_dataset and epoch % 5 == 0:
            print("\nRunning validation...", flush=True)
            model.eval()
            val_total, val_l1, val_hist, val_mass, val_spectral, val_hd = 0, 0, 0, 0, 0, 0
            
            with torch.no_grad():
                for val_batch in val_loader:
                    if len(val_batch) == 3:
                        val_x, val_y, _ = val_batch
                        if args.in_channel == 1:
                            val_x = val_x.unsqueeze(2)
                            val_y = val_y.unsqueeze(1)
                    else:
                        val_x, val_y = val_batch
                    
                    val_x, val_y = val_x.to(DEVICE), val_y.to(DEVICE)
                    val_delta = val_y-val_x
                    if use_amp:
                        with torch.amp.autocast(device_type=str(DEVICE)):
                            val_pred = model(val_x)
                        pred = val_pred.float()
                        val_delta_f = val_delta.float()
                    else:
                        pred = model(val_x)
                        val_delta_f = val_delta
                    
                    l1, hist_l, spectral_l, hd_l, mass_l = return_total_loss_multichannel(val_delta_f, pred, config, input_state=val_x)

                    losses = torch.stack([l1, hist_l, mass_l, spectral_l, hd_l])
                    w_t = torch.tensor(w_print, device=DEVICE)
                    total_loss = (w_t * losses).sum()
                    
                    batch_size = val_x.size(0)
                    val_total += total_loss.item() * batch_size
                    val_l1 += l1.item() * batch_size
                    val_hist += hist_l.item() * batch_size
                    val_mass += mass_l.item() * batch_size
                    val_spectral += spectral_l.item() * batch_size
                    val_hd += hd_l.item() * batch_size
            
            N_val = len(val_dataset)
            val_loss_history.append([
                val_total / N_val, val_l1 / N_val, val_hist / N_val,
                val_mass / N_val, val_spectral / N_val, val_hd / N_val
            ])
            stop_early = Earlystopper(val_total / N_val)
            
            print(f"Validation: Total={val_total/N_val:.4f} | "
                  f"L1={val_l1/N_val:.4f} HD={val_hd/N_val:.4f} "
                  f"Mass={val_mass/N_val:.4f} Hist={val_hist/N_val:.4f} "
                  f"Spectral={val_spectral/N_val:.4f}\n", flush=True)
            
            # TensorBoard validation logging
            writer.add_scalar("Val_Loss/Total", val_total / N_val, epoch)
            writer.add_scalar("Val_Loss/L1", val_l1 / N_val, epoch)
            writer.add_scalar("Val_Loss/Hist", val_hist / N_val, epoch)
            writer.add_scalar("Val_Loss/Mass", val_mass / N_val, epoch)
            writer.add_scalar("Val_Loss/Spectral", val_spectral / N_val, epoch)
            writer.add_scalar("Val_Loss/HighDensity", val_hd / N_val, epoch)
        
        # Periodic checkpointing at phase boundaries
        checkpoint_epochs = [24, 49, 69, 79, 89]  # End of each phase
        if epoch in checkpoint_epochs:
            print(f"\nSaving checkpoint at epoch {epoch+1}...", flush=True)
            torch.save({
                'epoch': epoch + 1,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'loss': loss_history,
                'val_loss': val_loss_history
            }, checkpoint_path)
        
        if stop_early:
            print(f"Early stopping triggered at epoch {epoch+1}.", flush=True)
            break
                
    writer.close()
    loss_history = np.array(loss_history)
    if dont_save:
        return min(loss_history[:, 0])  # Return minimum loss
    # Final save
    print("\nSaving final model...", flush=True)
    torch.save({
        'epoch': args.epochs,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'scheduler_state_dict': scheduler.state_dict(),
        'loss': loss_history,
        'val_loss': val_loss_history
    }, checkpoint_path)
    
    # Plot training curves
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.plot(loss_history[:, 0], marker='o', label='Total Loss')
    plt.plot(loss_history[:, 1], marker='.', label='L1 Loss')
    plt.plot(loss_history[:, 2], marker='*', label='Histogram Loss')
    plt.plot(loss_history[:, 3], marker='+', label='Mass Loss')
    plt.plot(loss_history[:, 4], marker='x', label='Spectral Loss')
    plt.plot(loss_history[:, 5], marker='1', label='High Density Loss')
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.xlabel("Epoch")
    plt.ylabel("Avg Weighted Loss")
    plt.title(f"{args.run_name} Training Loss")
    plt.savefig(f"models/plots/training_loss_{args.run_name}.png", dpi=150)
    plt.close()
    
    # Plot validation curves if available
    if val_loss_history:
        val_loss_history = np.array(val_loss_history)
        plt.figure(figsize=(10, 6))
        plt.plot(val_loss_history[:, 0], marker='o', label='Total Loss')
        plt.plot(val_loss_history[:, 1], marker='.', label='L1 Loss')
        plt.plot(val_loss_history[:, 2], marker='*', label='Histogram Loss')
        plt.plot(val_loss_history[:, 3], marker='+', label='Mass Loss')
        plt.plot(val_loss_history[:, 4], marker='x', label='Spectral Loss')
        plt.plot(val_loss_history[:, 5], marker='1', label='High Density Loss')
        plt.yscale('log')
        plt.legend()
        plt.grid()
        plt.tight_layout()
        plt.xlabel("Validation Check (every 5 epochs)")
        plt.ylabel("Avg Weighted Loss")
        plt.title(f"{args.run_name} Validation Loss")
        plt.savefig(f"models/plots/validation_loss_{args.run_name}.png", dpi=150)
        plt.close()
    
    total_time = time.perf_counter() - training_start
    print(f"\nTraining completed in {str(datetime.timedelta(seconds=int(total_time)))}", flush=True)
    
    return min(loss_history[:, 0])  # Return minimum loss
