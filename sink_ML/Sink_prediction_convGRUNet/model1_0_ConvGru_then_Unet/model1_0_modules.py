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
from torch.utils.data import DataLoader, TensorDataset
import json
from torchvision.utils import save_image
import logging
from torch.utils.tensorboard import SummaryWriter

logging.basicConfig(format="%(asctime)s - %(levelname)s: %(message)s", level=logging.INFO, datefmt="%I:%M:%S")

DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
#DEVICE = torch.device('cpu')

MODELFILE = 'nnmodel.dict'


def count_parameters(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)

def print_memory_stats():
    if torch.cuda.is_available():
        print(f"Allocated: {torch.cuda.memory_allocated()/1024**3:.2f} GB", flush=True)
        print(f"Reserved: {torch.cuda.memory_reserved()/1024**3:.2f} GB", flush=True)

import os
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

import time, datetime, torch
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

    est_remain_epoch_s = batch_time_ema * batches_left

    # compute average epoch time so far (approx)
    completed_epochs = max(1, (epoch - start_epoch + 1))
    elapsed_since_training_start = time.perf_counter() - training_start
    avg_epoch_time = elapsed_since_training_start / completed_epochs
    epochs_left = max(0, args.epochs - epoch - 1)
    est_total_remaining_s = est_remain_epoch_s + epochs_left * avg_epoch_time

    def fmt_s(sec):
        return str(datetime.timedelta(seconds=int(sec)))

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
    if (batches_done % 500) == 0:
        print()

    return batch_time_ema

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

def get_all_images(input_dir):
    output_dic = {}
    for filename in sorted(os.listdir(input_dir)):
        infile = open(os.path.join(input_dir, filename), 'r')
        i_file = json.load(infile)
        num = filename.split('_')[1]
        img_dic = i_file['img']
        label = i_file['label']
        output_dic[num] = {'img': img_dic, 'label': label}
        infile.close()
    return output_dic

def get_subset_images(input_dir, indices):
    output_dic = {}
    all_files = sorted(os.listdir(input_dir))
    for i in indices:
        filename = all_files[i]
        infile = open(os.path.join(input_dir, filename), 'r')
        i_file = json.load(infile)
        num = filename.split('_')[1]
        img_dic = i_file['img']
        label = i_file['label']
        output_dic[num] = {'img': img_dic, 'label': label}
        infile.close()
    return output_dic

from torch.utils.data import Dataset, DataLoader
import h5py
# Custom Dataset class for HDF5 data
class HDF5Dataset(Dataset):
    """
    Custom PyTorch Dataset for loading data from HDF5 files efficiently.
    This allows for lazy loading - data is only loaded when needed.
    """
    def __init__(self, hdf5_file, indices=None):
        """
        Initialize HDF5Dataset
        
        Parameters:
        -----------
        hdf5_file : str
            Path to HDF5 file
        indices : list, optional
            Specific indices to use. If None, uses all available data.
        """
        self.hdf5_file = hdf5_file
        
        # Get available indices from the file
        with h5py.File(hdf5_file, 'r') as f:
            all_keys = list(f['input_images'].keys())
            
        if indices is not None:
            # Convert indices to corresponding keys
            self.keys = [all_keys[i] for i in indices if i < len(all_keys)]
        else:
            self.keys = all_keys
            
    def __len__(self):
        return len(self.keys)
    
    def __getitem__(self, idx):
        key = self.keys[idx]
        
        with h5py.File(self.hdf5_file, 'r') as f:
            # Load input and output images
            input_img = torch.tensor(f['input_images'][key][...], dtype=torch.float32)
            output_img = torch.tensor(f['output_images'][key][...], dtype=torch.float32)
            
            # Load label
            label_str = f['labels'][key][()]
            label_str = label_str.decode('utf-8')
            label_str = eval(label_str)  # Convert string representation of list to actual list
        return input_img, output_img, label_str

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


class AttentionGate(nn.Module):
    def __init__(self, F_g, F_l, F_int, args):
        super().__init__()
        # Reduce gating signal channels
        self.W_g = nn.Sequential(
            nn.Conv3d(F_g, F_int, kernel_size=1, bias=True, padding_mode='circular' ),
            nn.BatchNorm3d(F_int)
        )
        # Reduce encoder feature channels
        self.W_x = nn.Sequential(
            nn.Conv3d(F_l, F_int, kernel_size=1, bias=True, padding_mode='circular'),
            nn.BatchNorm3d(F_int)
        )
        # Combine & compute attention coefficients
        self.psi = nn.Sequential(
            nn.Conv3d(F_int, 1, kernel_size=1, bias=True, padding_mode='circular'),
            nn.BatchNorm3d(1),
            nn.Sigmoid()
        )
        self.relu = nn.LeakyReLU(negative_slope=args.LRLU_slope, inplace=True)

    def forward(self, g, x):
        """
        x: encoder feature map (F_l channels)
        g: gating signal    (F_g channels) from decoder
        """
        g1 = self.W_g(g)             # [B, F_int, D, H, W]
        x1 = self.W_x(x)             # [B, F_int, D, H, W]
        psi = self.relu(g1 + x1)     # broadcast-sum
        psi = self.psi(psi)          # [B, 1, D, H, W]
        return x * psi               # gated encoder features


class Convolution_Block(nn.Module):
    def __init__(self,ch_in,ch_out, args):
        super().__init__()
        self.conv = nn.Sequential(
            nn.Conv3d(ch_in, ch_out, kernel_size=args.conv_param[0],stride=args.conv_param[1],padding=args.conv_param[2],
                       dilation=args.conv_param[3], bias=True, padding_mode='circular'),
            nn.BatchNorm3d(ch_out),
            nn.LeakyReLU(negative_slope=args.LRLU_slope, inplace=True),
            nn.Dropout3d(p=args.dropout_rate),
            nn.Conv3d(ch_out, ch_out, kernel_size=args.conv_param[0],stride=args.conv_param[1],padding=args.conv_param[2],
                      dilation=args.conv_param[3], bias=True, padding_mode='circular'),
            nn.BatchNorm3d(ch_out),
            nn.LeakyReLU(negative_slope=args.LRLU_slope, inplace=True),
            nn.Dropout3d(p=args.dropout_rate)
        )

    def forward(self,x):
        x = self.conv(x)
        return x

class Upsampling_Block(nn.Module):
    def __init__(self,ch_in,ch_out, args):
        super().__init__()
        self.up = nn.Sequential(
            nn.Upsample(scale_factor=2),
            nn.Conv3d(ch_in,ch_out,kernel_size=args.conv_param[0],stride=args.conv_param[1],padding=args.conv_param[2],
                      dilation=args.conv_param[3], bias=True, padding_mode='circular'),
		    nn.BatchNorm3d(ch_out),
			nn.LeakyReLU(negative_slope=args.LRLU_slope, inplace=True)
        )

    def forward(self,x):
        x = self.up(x)
        return x


class UNet3D(nn.Module):
    def __init__(self, in_ch, out_ch, args):
        """
        A 3D U-Net with attention gates, 4 level deep
        Args:
            in_ch: number of input channels (e.g., 1 for grayscale)
            out_ch: number of output channels (e.g., 1 for grayscale)
            base_f: base number of filters in the first layer
        """
        super().__init__()

        self.Maxpool = nn.MaxPool3d(kernel_size=2,stride=2)

        self.Conv1 = Convolution_Block(ch_in=in_ch,ch_out=args.base_f, args=args)
        self.Conv2 = Convolution_Block(ch_in=args.base_f,ch_out=args.base_f*2, args=args)
        if args.depth >= 3:
            self.Conv3 = Convolution_Block(ch_in=args.base_f*2,ch_out=args.base_f*4, args=args)
        if args.depth >= 4:
            self.Conv4 = Convolution_Block(ch_in=args.base_f*4,ch_out=args.base_f*8, args=args)
        if args.depth >= 5:
            self.Conv5 = Convolution_Block(ch_in=args.base_f*8,ch_out=args.base_f*16, args=args)

        # Upsampling blocks with attention gates
        if args.depth >= 5:
            self.Up4 = Upsampling_Block(ch_in=args.base_f*16,ch_out=args.base_f*8, args=args)
            self.Att4 = AttentionGate(F_g=args.base_f*8,F_l=args.base_f*8,F_int=args.base_f*4//args.attention_int_division, args=args)
            self.Up_conv4 = Convolution_Block(ch_in=args.base_f*16, ch_out=args.base_f*8, args=args)
        if args.depth >= 4:
            self.Up3 = Upsampling_Block(ch_in=args.base_f*8,ch_out=args.base_f*4, args=args)
            self.Att3 = AttentionGate(F_g=args.base_f*4,F_l=args.base_f*4,F_int=(args.base_f*4)//args.attention_int_division, args=args)
            self.Up_conv3 = Convolution_Block(ch_in=args.base_f*8, ch_out=args.base_f*4, args=args)
        if args.depth >= 3:
            self.Up2 = Upsampling_Block(ch_in=args.base_f*4,ch_out=args.base_f*2, args=args)
            self.Att2 = AttentionGate(F_g=args.base_f*2,F_l=args.base_f*2,F_int=(args.base_f*2)//args.attention_int_division, args=args)
            self.Up_conv2 = Convolution_Block(ch_in=args.base_f*4, ch_out=args.base_f*2, args=args)

        self.Up1 = Upsampling_Block(ch_in=args.base_f*2,ch_out=args.base_f, args=args)
        self.Att1 = AttentionGate(F_g=args.base_f,F_l=args.base_f,F_int=args.base_f//args.attention_int_division, args=args)
        self.Up_conv1 = Convolution_Block(ch_in=args.base_f*2, ch_out=args.base_f, args=args)

        self.Conv_1x1 = nn.Conv3d(args.base_f, out_ch, kernel_size=1, stride=1, padding=0)

    def forward(self,x):
        # encoding path
        x0 = self.Conv1(x)

        #First Encoder
        x1 = self.Maxpool(x0)
        x1 = self.Conv2(x1)

        if hasattr(self, 'Conv3'):
            #Second Encoder
            x2 = self.Maxpool(x1)
            x2 = self.Conv3(x2)
        if hasattr(self, 'Conv4'):
            #Third Encoder
            x3 = self.Maxpool(x2)
            x3 = self.Conv4(x3)
        if hasattr(self, 'Conv5'):
            #Fourth Encoder
            x4 = self.Maxpool(x3)
            x4 = self.Conv5(x4)

        #decoding + concat path
        if hasattr(self, 'Conv5'):
            d3 = self.Up4(x4)
            a4 = self.Att4(g=d3,x=x3)
            d3 = torch.cat((a4,d3),dim=1)
            x3 = self.Up_conv4(d3)

        if hasattr(self, 'Conv4'):
            d2 = self.Up3(x3)
            a3 = self.Att3(g=d2,x=x2)
            d2 = torch.cat((a3,d2),dim=1)
            x2 = self.Up_conv3(d2)

        if hasattr(self, 'Conv3'):
            d1 = self.Up2(x2)
            a2 = self.Att2(g=d1,x=x1)
            d1 = torch.cat((a2,d1),dim=1)
            x1 = self.Up_conv2(d1)

        d0 = self.Up1(x1)
        a1 = self.Att1(g=d0,x=x0)
        d0 = torch.cat((a1,d0),dim=1)
        d0 = self.Up_conv1(d0)

        # Final 1x1 convolution to get the output channels
        return self.Conv_1x1(d0)


class hybridConvGRU3DNET(nn.Module):
    def __init__(self, args, argsGRU):
        """
        Initialize the hybridConvGRU3DNET model
        :param args:
            args for the UNet3D model
        :param in_seq: int
            Number of input sequences, how many time states given
        :param out_seq: int
            Number of output sequences, how many time states to predict
        :param num_layers: int
            number of layers in GRU
        :param kernels: (int,int,...)
            kernel size for each layer
        :param hidden_dims: (int,int,...)
            channel size for each layer
        """
        super().__init__()

        self.ConvGRUNet = ConvGRUNet3D(input_channels=argsGRU.in_channel, output_channels =argsGRU.in_channel, image_size=args.image_size, in_seq = argsGRU.in_seq, out_seq = argsGRU.out_seq, 
                                       num_layers = argsGRU.num_layers, kernels = argsGRU.kernels, hidden_dims = argsGRU.hidden_dims, bias = argsGRU.bias)
        self.UNet3D = UNet3D(in_ch=argsGRU.in_channel, out_ch=argsGRU.out_channel, args=args)

    def forward(self, input_tensor):
        """
        :param input_tensor: (b, t_in, c, l)
        :param hidden_state:
        :return: output_tensor: (b, t_out, c, l)
        """
        GRU_output = self.ConvGRUNet(input_tensor)  #(b, t_out, c, l, h ,w)
        unet_output = self.UNet3D(GRU_output[:,0,:,:,:,:]) #(b, c, l, h, w)
        return unet_output
        

class ConvGRUNet3D(nn.Module):
    def __init__(self, input_channels=1, output_channels = 1, image_size=64, in_seq = 2, out_seq = 1, num_layers = 1, kernels = [3], hidden_dims = [16], bias = True):
        """
        Initialize the ConvGRU model
        :param input_channels: int
            Number of input channels
        :param in_seq: int
            Number of input sequences, how many time states given
        :param out_seq: int
            Number of output sequences, how many time states to predict
        :param num_layers: int
            number of layers in GRU
        :param kernels: (int,int,...)
            kernel size for each layer
        :param hidden_dims: (int,int,...)
            channel size for each layer
        :param bias: bool
            Whether or not to add the bias.
        """
        super().__init__()
        self.input_channels = input_channels
        self.output_channels = output_channels
        self.image_size = image_size
        self.in_seq = in_seq
        self.out_seq = out_seq
        self.num_layers = num_layers
        self.kernels = kernels
        self.hidden_dims = hidden_dims
        self.bias = bias


        cell_list = []
        #check condition
        if not (len(self.kernels)==len(self.hidden_dims) and len(self.kernels)==num_layers):
            raise ValueError("kernel and hidden layer lists don't match number of layers!")
        
        for i in range(0, self.num_layers):
            cur_input_dim = self.input_channels if i == 0 else self.hidden_dims[i - 1]
            cell_list.append(ConvGRUCell3D(input_channels=cur_input_dim,
                                         image_size=self.image_size,
                                         hidden_dim=self.hidden_dims[i],
                                         kernel_size=self.kernels[i],
                                         bias=self.bias))

        # convert python list to pytorch module
        self.cell_list = nn.ModuleList(cell_list)

        self.out_cnn = nn.Conv3d(in_channels=self.hidden_dims[-1],
                                    out_channels=self.output_channels,
                                    kernel_size=self.kernels[-1],
                                    padding=self.kernels[-1]//2,
                                    bias=self.bias,
                                    padding_mode='circular')

        self.mse=nn.MSELoss()
        self.log_derivative_weight = nn.Parameter(torch.tensor(0.0)) 
        self.hl = nn.HuberLoss(delta=0.2)
        self.l1 = nn.L1Loss()

    def criterion(self,guess,target, initial=None):
        L1 = self.l1(target,guess)
        return L1

    def forward(self, input_tensor, hidden_state=None):
        """
        :param input_tensor: (b, t_in, c, l, h, w)
        :param hidden_state:
        :return: output_tensor: (b, t_out, c, l, h, w)
        """
        hidden_state = self._init_hidden(batch_size=input_tensor.size(0))
        device = input_tensor.device
        hidden_state = [h.to(device) for h in hidden_state]

        layer_output_list = []
        cur_layer_input = input_tensor
        for layer_idx in range(self.num_layers):
            h = hidden_state[layer_idx]
            output_inner = []
            for t in range(self.in_seq):
                # input current hidden and cell state then compute the next hidden and cell state through ConvGRU forward function
                h = self.cell_list[layer_idx](input_tensor=cur_layer_input[:, t, :, :, :, :], #(b,t,c,l, h, w)
                                              h_cur=h)
                output_inner.append(h)

            layer_output = torch.stack(output_inner, dim=1)
            cur_layer_input = layer_output

            layer_output_list.append(layer_output)
        
        last_layer_output = layer_output_list[-1] # shape (b,t_in,hi,l,h,w)
        out_arr = last_layer_output[:,-1,:,:,:,:] # shape (b,h,l)
        out_arr = self.out_cnn(out_arr)[:,None,:,:,:,:] # shape (b,t_out,c,l)

        return out_arr

    def _init_hidden(self, batch_size):
        init_states = []
        for i in range(self.num_layers):
            init_states.append(self.cell_list[i].init_hidden(batch_size))
        return init_states


class ConvGRUCell3D(nn.Module):
    def __init__(self, input_channels, image_size, hidden_dim, kernel_size, bias):
        """
        Initialize a ConvGRU cell
        :param input_channels: int
            number of input channels
        :param hidden_dim: int
            Number of channels of hidden state
        :param kernel_size: (int, int)
            Size of the convolutional kernel.
        :param bias: bool
            Whether or not to add the bias.
        """
        super(ConvGRUCell3D, self).__init__()
        self.input_channels = input_channels
        self.image_size = image_size
        self.padding = kernel_size // 2
        self.hidden_dim = hidden_dim
        self.bias = bias

        self.conv_gates = nn.Conv3d(in_channels=self.input_channels + self.hidden_dim,
                                    out_channels=2*self.hidden_dim,  # for update_gate,reset_gate respectively
                                    kernel_size=kernel_size,
                                    padding=self.padding,
                                    bias=self.bias,
                                    padding_mode='circular')

        self.conv_can = nn.Conv3d(in_channels=self.input_channels + self.hidden_dim,
                              out_channels=self.hidden_dim, # for candidate neural memory
                              kernel_size=kernel_size,
                              padding=self.padding,
                              bias=self.bias,
                              padding_mode='circular')

    def init_hidden(self, batch_size):
        return torch.zeros(batch_size, self.hidden_dim, self.image_size, self.image_size, self.image_size)

    def forward(self, input_tensor, h_cur):
        """
        :param self:
        :param input_tensor: (b, c, l)
            input is actually the target_model
        :param h_cur: (b, c_hidden, l)
            current hidden state
        :return: h_next
            next hidden state
        """
        combined = torch.cat([input_tensor, h_cur], dim=1)
        combined_conv = self.conv_gates(combined)

        gamma, beta = torch.split(combined_conv, self.hidden_dim, dim=1)
        reset_gate = torch.sigmoid(gamma)
        update_gate = torch.sigmoid(beta)

        combined = torch.cat([input_tensor, reset_gate*h_cur], dim=1)
        cc_cnm = self.conv_can(combined)
        cnm = torch.tanh(cc_cnm)

        h_next = (1 - update_gate) * h_cur + update_gate * cnm
        return h_next



def compute_spectral_loss(pred, target):
    #pred, target are shaped [B, 1, D, H, W]
    pred_fft = torch.fft.fftn(pred.squeeze(1), dim=(1,2,3))   # [B, D, H, W], complex
    targ_fft = torch.fft.fftn(target.squeeze(1), dim=(1,2,3))
    #find the absolute magnitude, whihc gives the power spectrum
    pred_mag = torch.abs(pred_fft)
    targ_mag = torch.abs(targ_fft)
    log_pred = torch.log1p(pred_mag)
    log_true = torch.log1p(targ_mag)
    loss = F.mse_loss(log_pred, log_true)
    return loss

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

from torch.utils.tensorboard import SummaryWriter
#FIX DELTA ISSUES
def train_unet_delta(input_arr, output_arr, labels, args, argsGRU):
    X = torch.stack(input_arr, dim=0).float().unsqueeze(2)  # [N,t, 1,64, 64, 64]
    Y = torch.stack(output_arr, dim=0).float().unsqueeze(1) # [N, 1,64, 64, 64]

    dataset = TensorDataset(X, Y)
    loader  = DataLoader(dataset,
                         batch_size=args.batch_size,
                         shuffle=True,
                         num_workers=2,    # adjust based on your CPU
                         pin_memory=False)  # if using GPU

    # Initialize model, loss function, and optimizer
    model = hybridConvGRU3DNET(args, argsGRU).to(DEVICE, non_blocking=True)
    # For trainable weights
    n_losses = 5  # L1, hist, Mass, Spectral, High Density
    loss_w = torch.nn.Parameter(torch.ones(n_losses, device=DEVICE), requires_grad=True)
    optimizer = optim.AdamW(list(model.parameters()) + [loss_w], lr=args.lr, weight_decay=args.Adamw_weight_decay)
    writer = SummaryWriter(log_dir=f"models/logs/{args.run_name}")

    start_epoch = 0
    loss_history = []
    #check if model already exists
    if os.path.exists(f"models/{args.run_name}_{MODELFILE}"):
        checkpoint = torch.load(f"models/{args.run_name}_{MODELFILE}", map_location=DEVICE, weights_only=False)
        model.load_state_dict(checkpoint['model_state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        start_epoch = checkpoint['epoch']
        loss_history = checkpoint['loss']
        print(f"Model loaded for training with {count_parameters(model):,} trainable parameters.", flush=True)
    else:
        print(f"Model initialized for training with {count_parameters(model):,} trainable parameters.", flush=True)
    
    L1_scale = 5.0  # scale L1 loss to be more influential
    hist_scale = 10.0  # scale histogram loss
    mass_scale = 0.1   # give mass only 10% initial influence
    spectral_scale = 0.5  # scale spectral loss
    hd_scale = 50.0 #high density loss scale


    if args.loss_type == 'static':
        static_w = [1.0, 1.0, 1.0, 1.0, 1.0]  # initial static weights for L1, hist, Mass, Spectral, High Density
    if args.loss_type == 'dynamic':
        prev_losses = [1.0]*n_losses    # dummy for epoch 0
        prev2_losses = [1.0]*n_losses   # dummy for epoch -1
        w = [1.0]*n_losses

    # === Training loop ===
    for epoch in range(start_epoch, args.epochs):
        model.train()
        total_epoch_loss = 0.0
        sum_l1, sum_hist, sum_mass, sum_spectral, sum_hd = 0.0, 0.0, 0.0, 0.0, 0.0

        for x, y in loader:
            x, y = x.to(DEVICE, non_blocking=True), y.to(DEVICE, non_blocking=True)
            delta = y-x
            optimizer.zero_grad()
            deltapred = model(x)
            #Computing individual losses
            #L1 loss
            l1 = F.l1_loss(deltapred, delta)
            l1 = l1 * L1_scale  # scale L1 loss
            #Histogram Loss
            hist_l = hist_loss(deltapred, delta, bins=32)
            hist_l = hist_l * hist_scale  # scale histogram loss
            #Mass loss
            mass_t = delta.sum()
            mass_p = deltapred.sum()
            mass_error = torch.abs(mass_p - mass_t)
            mass_l = torch.log1p(mass_error)  # log1p to avoid large values
            mass_l = mass_l * mass_scale  # scale mass loss
            # Spectral loss using power spectrum
            spectral_l = compute_spectral_loss(deltapred, delta)
            spectral_l = spectral_l * spectral_scale  # scale spectral loss
            #High Density loss using heaviside function as filter
            with torch.no_grad():
                hd_true = torch.heaviside(delta-torch.quantile(delta.flatten(), 0.99), torch.tensor([1.0]))
                hd_pred = torch.heaviside(deltapred-torch.quantile(deltapred.flatten(), 0.99), torch.tensor([1.0]))
            hd_loss = F.l1_loss(hd_pred, hd_true)
            hd_l = hd_loss*hd_scale

            #Stacking losses
            losses = torch.stack([l1, hist_l, mass_l, spectral_l, hd_l])
            if args.loss_type == 'static':
                total_loss = (static_w * losses).sum()
            else:
                w_t = torch.tensor(w, device=DEVICE)
                total_loss = (w_t * losses).sum()

            total_loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)  # optional
            optimizer.step()
            # accumulate
            batch_size = x.size(0)
            total_epoch_loss += total_loss.item() * batch_size
            sum_l1     += l1.item()      * batch_size
            sum_hist    += hist_l.item()     * batch_size
            sum_mass   += mass_l.item()  * batch_size
            sum_spectral   += spectral_l.item()  * batch_size
            sum_hd   += hd_l.item()  * batch_size
            

        #compute per-sample averages
        N = len(dataset)
        avg_total = total_epoch_loss / N
        avg_l1    = sum_l1     / N
        avg_hist  = sum_hist   / N
        avg_mass  = sum_mass   / N
        avg_spectral = sum_spectral / N
        avg_hd = sum_hd / N
        if args.loss_type == 'dynamic':
            avg_losses = [avg_l1, avg_hist, avg_mass, avg_spectral, avg_hd]
            if epoch >= 2:
                #compute r_i = L_i(t-1)/L_i(t-2)
                r = [prev_losses[i]/prev2_losses[i] for i in range(n_losses)]
                #DWA weights
                T = args.DWA_temperature #defined as temperature in the paper, smoothes extremes
                K = n_losses   #defined as number of tasks in the paper
                exp_r = [np.exp(r_i/T) for r_i in r]
                w = [(K * e) / sum(exp_r) for e in exp_r]
            else:
                # use static equal weights for first two epochs
                w = [1.0]*n_losses

            prev2_losses = prev_losses
            prev_losses  = avg_losses
            w_print = w
        else:
            w_print = static_w #for console

        loss_history.append(float(avg_total))

        # log to console
        print(f"Epoch {epoch+1}: "
              f"Total={avg_total:.4f} | "
              f"L1={avg_l1:.4f} "
              f"HighDensity={avg_hd:.4f} "
              f"Mass={avg_mass:.4f}  Hist={avg_hist:.4f} "
              f"Spectral={avg_spectral:.4f} "
              f"Weights: {*w_print,}")

        # log to TensorBoard
        writer.add_scalar("Loss/Total",    avg_total, epoch)
        writer.add_scalar("Loss/L1",       avg_l1,    epoch)
        writer.add_scalar("Loss/Hist",     avg_hist,   epoch)
        writer.add_scalar("Loss/Mass",     avg_mass,  epoch)
        writer.add_scalar("Loss/Spectral", avg_spectral, epoch)
        writer.add_scalar("Loss/HighDensity", avg_hd, epoch)

        if args.loss_type!='static':
            # log the learned weights (after softmax)
            curr_w = torch.softmax(loss_w, dim=0).detach().cpu().tolist()
            for i, name in enumerate(["L1","Hist","Mass", "Spectral", "HighDensity"]):
                writer.add_scalar(f"Weights/{name}", curr_w[i], epoch)

    writer.close()

    # === Save & plot ===
    torch.save({'epoch': args.epochs, 'model_state_dict': model.state_dict(), 'optimizer_state_dict': optimizer.state_dict(), 'loss': loss_history},
                f"models/{args.run_name}_{MODELFILE}")
    if args.loss_type != 'static':
        torch.save(loss_w.detach().cpu(), f"models/plots/loss_data/{args.run_name}_loss.pt")
    import matplotlib.pyplot as plt
    print(loss_history)
    plt.plot(loss_history, marker='o')
    plt.xlabel("Epoch")
    plt.ylabel("Avg Weighted Loss")
    plt.title(f"{args.run_name} Training Loss")
    plt.savefig(f"models/plots/training_loss/{args.run_name}.png")
    plt.close()

    return min(loss_history)  # return the minimum loss for this run

def train_unet(input_arr, output_arr, labels, args, argsGRU):
    print(f"Training convGRUNet with {len(input_arr)} samples, {args.run_name} run", flush=True)
    X = torch.stack(input_arr, dim=0).float().unsqueeze(2)  # [N, t, 1,64, 64, 64]
    Y = torch.stack(output_arr, dim=0).float().unsqueeze(1) # [N, 1,64, 64, 64]

    num_workers = args.num_workers if hasattr(args, 'num_workers') else min(8, os.cpu_count()//2)
    print(f"Using {num_workers} DataLoader workers", flush=True)
    dataset = TensorDataset(X, Y)
    loader  = DataLoader(dataset,
                         batch_size=args.batch_size,
                         shuffle=True,
                         num_workers=num_workers,    # adjust based on your CPU
                         pin_memory=True,  # if using GPU
                        persistent_workers=(num_workers>0))
    # Initialize model, loss function, and optimizer
    model = hybridConvGRU3DNET(args, argsGRU).to(DEVICE, non_blocking=True)
    # For trainable weights
    n_losses = 5  # L1, hist, Mass, Spectral, High Density
    loss_w = torch.nn.Parameter(torch.ones(n_losses, device=DEVICE), requires_grad=True)
    optimizer = optim.AdamW(list(model.parameters()) + [loss_w], lr=args.lr, weight_decay=args.Adamw_weight_decay)
    writer = SummaryWriter(log_dir=f"models/logs/{args.run_name}")

    start_epoch = 0
    loss_history = []
    #check if model already exists
    if os.path.exists(f"models/{args.run_name}_{MODELFILE}"):
        checkpoint = torch.load(f"models/{args.run_name}_{MODELFILE}", map_location=DEVICE, weights_only=False)
        model.load_state_dict(checkpoint['model_state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        start_epoch = checkpoint['epoch']
        loss_history = checkpoint['loss']
        print(f"Model loaded for training with {count_parameters(model):,} trainable parameters.", flush=True)
    else:
        print(f"Model initialized for training with {count_parameters(model):,} trainable parameters.", flush=True)
    
    L1_scale = 5.0  # scale L1 loss to be more influential
    hist_scale = 1.0  # scale histogram loss
    mass_scale = 0.1   # give mass only 10% initial influence
    spectral_scale = 0.5  # scale spectral loss
    hd_scale = 1.0 #high density loss scale


    if args.loss_type == 'static':
        static_w = [1.0, 1.0, 1.0, 1.0, 1.0]  # initial static weights for L1, hist, Mass, Spectral, High Density
    if args.loss_type == 'dynamic':
        prev_losses = [1.0]*n_losses    # dummy for epoch 0
        prev2_losses = [1.0]*n_losses   # dummy for epoch -1
        w = [1.0]*n_losses

    # === Training loop ===
    for epoch in range(start_epoch, args.epochs):
        model.train()
        total_epoch_loss = 0.0
        sum_l1, sum_hist, sum_mass, sum_spectral, sum_hd = 0.0, 0.0, 0.0, 0.0, 0.0

        for x, y in loader:
            x, y = x.to(DEVICE, non_blocking=True), y.to(DEVICE, non_blocking=True)
            optimizer.zero_grad()
            pred = model(x)
            #Computing individual losses
            #L1 loss
            l1 = F.l1_loss(y, pred)
            l1 = l1 * L1_scale  # scale L1 loss
            #Histogram Loss
            hist_l = hist_loss(y, pred, bins=32)
            hist_l = hist_l * hist_scale  # scale histogram loss
            #Mass loss
            mass_t = torch.pow(10, y).sum()
            mass_p = torch.pow(10, pred).sum()
            mass_error = torch.abs(mass_p - mass_t) / (mass_t + 1e-8)
            mass_l = torch.log1p(mass_error)  # log1p to avoid large values
            mass_l = mass_l * mass_scale  # scale mass loss
            # Spectral loss using power spectrum
            spectral_l = compute_spectral_loss(y, pred)
            spectral_l = spectral_l * spectral_scale  # scale spectral loss
            #High Density loss using heaviside function as filter
            with torch.no_grad():
                hd_true = torch.heaviside(y-torch.quantile(y.flatten(), 0.99), torch.tensor([1.0]))
                hd_pred = torch.heaviside(pred-torch.quantile(pred.flatten(), 0.99), torch.tensor([1.0]))
            hd_loss = F.l1_loss(hd_pred, hd_true)
            hd_l = hd_loss*hd_scale

            #Stacking losses
            losses = torch.stack([l1, hist_l, mass_l, spectral_l, hd_l])
            if args.loss_type == 'static':
                total_loss = (static_w * losses).sum()
            else:
                w_t = torch.tensor(w, device=DEVICE)
                total_loss = (w_t * losses).sum()

            total_loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)  # optional
            optimizer.step()
            # accumulate
            batch_size = x.size(0)
            total_epoch_loss += total_loss.item() * batch_size
            sum_l1     += l1.item()      * batch_size
            sum_hist    += hist_l.item()     * batch_size
            sum_mass   += mass_l.item()  * batch_size
            sum_spectral   += spectral_l.item()  * batch_size
            sum_hd   += hd_l.item()  * batch_size
            

        #compute per-sample averages
        N = len(dataset)
        avg_total = total_epoch_loss / N
        avg_l1    = sum_l1     / N
        avg_hist  = sum_hist   / N
        avg_mass  = sum_mass   / N
        avg_spectral = sum_spectral / N
        avg_hd = sum_hd / N
        if args.loss_type == 'dynamic':
            avg_losses = [avg_l1, avg_hist, avg_mass, avg_spectral, avg_hd]
            if epoch >= 2:
                #compute r_i = L_i(t-1)/L_i(t-2)
                r = [prev_losses[i]/prev2_losses[i] for i in range(n_losses)]
                #DWA weights
                T = args.DWA_temperature #defined as temperature in the paper, smoothes extremes
                K = n_losses   #defined as number of tasks in the paper
                exp_r = [np.exp(r_i/T) for r_i in r]
                w = [(K * e) / sum(exp_r) for e in exp_r]
            else:
                # use static equal weights for first two epochs
                w = [1.0]*n_losses

            prev2_losses = prev_losses
            prev_losses  = avg_losses
            w_print = w
        else:
            w_print = static_w #for console

        loss_history.append([float(avg_total), float(avg_l1), float(avg_hist), float(avg_mass), float(avg_spectral), float(avg_hd)])

        # log to console
        print(f"Epoch {epoch+1}: "
              f"Total={avg_total:.4f} | "
              f"L1={avg_l1:.4f} "
              f"HighDensity={avg_hd:.4f} "
              f"Mass={avg_mass:.4f}  Hist={avg_hist:.4f} "
              f"Spectral={avg_spectral:.4f} "
              f"Weights: {*w_print,}")

        # log to TensorBoard
        writer.add_scalar("Loss/Total",    avg_total, epoch)
        writer.add_scalar("Loss/L1",       avg_l1,    epoch)
        writer.add_scalar("Loss/Hist",     avg_hist,   epoch)
        writer.add_scalar("Loss/Mass",     avg_mass,  epoch)
        writer.add_scalar("Loss/Spectral", avg_spectral, epoch)
        writer.add_scalar("Loss/HighDensity", avg_hd, epoch)

        if args.loss_type!='static':
            # log the learned weights (after softmax)
            curr_w = torch.softmax(loss_w, dim=0).detach().cpu().tolist()
            for i, name in enumerate(["L1","Hist","Mass", "Spectral", "HighDensity"]):
                writer.add_scalar(f"Weights/{name}", curr_w[i], epoch)

    writer.close()

    # === Save & plot ===
    torch.save({'epoch': args.epochs, 'model_state_dict': model.state_dict(), 'optimizer_state_dict': optimizer.state_dict(), 'loss': loss_history},
                f"models/{args.run_name}_{MODELFILE}")
    if args.loss_type != 'static':
        torch.save(loss_w.detach().cpu(), f"models/plots/loss_data_{args.run_name}_loss.pt")
    import matplotlib.pyplot as plt
    loss_history = np.array(loss_history)
    plt.plot(loss_history[:,0], marker='o', label='Total Loss')
    plt.plot(loss_history[:,1], marker='.', label='L1 Loss')
    plt.plot(loss_history[:,2], marker='*', label='Histogram Loss')
    plt.plot(loss_history[:,3], marker='+', label='Mass Loss')
    plt.plot(loss_history[:,4], marker='x', label='Spectral Loss')
    plt.plot(loss_history[:,5], marker='1', label='High Density Loss')
    plt.yscale('log')  # log scale for better visibility of loss changes
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.xlabel("Epoch")
    plt.ylabel("Avg Weighted Loss")
    plt.title(f"{args.run_name} Training Loss")
    plt.savefig(f"models/plots/training_loss_{args.run_name}.png")
    plt.close()

    return min(loss_history[:,0])  # return the minimum loss for this run


# Modified training function that can handle HDF5Dataset
def train_unet_with_dataset(dataset, args, argsGRU):
    torch.cuda.empty_cache()
    torch.backends.cudnn.benchmark = False
    print(f"Training convGRUNet with {len(dataset)} samples, {args.run_name} run", flush=True)
    
    num_workers = args.num_workers if hasattr(args, 'num_workers') else min(8, os.cpu_count()//2)
    print(f"Using {num_workers} DataLoader workers", flush=True)
       
    loader = DataLoader(dataset,
                       batch_size=args.batch_size,
                       shuffle=True,
                       num_workers=num_workers,
                       pin_memory=True,
                       persistent_workers=(num_workers>0))
    
    # Initialize model, loss function, and optimizer
    model = hybridConvGRU3DNET(args, argsGRU).to(DEVICE, non_blocking=True)
    
    # For trainable weights
    n_losses = 5  # L1, hist, Mass, Spectral, High Density
    loss_w = torch.nn.Parameter(torch.ones(n_losses, device=DEVICE), requires_grad=True)
    optimizer = optim.AdamW(list(model.parameters()) + [loss_w], lr=args.lr, weight_decay=args.Adamw_weight_decay)
    writer = SummaryWriter(log_dir=f"models/logs/{args.run_name}")

    start_epoch = 0
    loss_history = []
    
    # Check if model already exists
    if os.path.exists(f"models/{args.run_name}_{MODELFILE}"):
        checkpoint = torch.load(f"models/{args.run_name}_{MODELFILE}", map_location=DEVICE, weights_only=False)
        model.load_state_dict(checkpoint['model_state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        start_epoch = checkpoint['epoch']
        loss_history = checkpoint['loss']
        print(f"Model loaded for training with {count_parameters(model):,} trainable parameters.", flush=True)
    else:
        print(f"Model initialized for training with {count_parameters(model):,} trainable parameters.", flush=True)
    
    # Loss scales
    L1_scale = 5.0
    hist_scale = 1.0
    mass_scale = 0.1
    spectral_scale = 0.5
    hd_scale = 1.0

    if args.loss_type == 'static':
        static_w = [1.0, 1.0, 1.0, 1.0, 1.0]
    if args.loss_type == 'dynamic':
        prev_losses = [1.0]*n_losses
        prev2_losses = [1.0]*n_losses
        w = [1.0]*n_losses

    # Training loop
    import time, datetime
    training_start = time.perf_counter()
    batch_time_ema = None
    ema_alpha = 0.15
    for epoch in range(start_epoch, args.epochs):
        model.train()
        total_epoch_loss = 0.0
        sum_l1, sum_hist, sum_mass, sum_spectral, sum_hd = 0.0, 0.0, 0.0, 0.0, 0.0
        torch.cuda.empty_cache()
        for batch_idx, batch_data in enumerate(loader):
            batch_start = time.perf_counter()
            if len(batch_data) == 3:
                x, y, labels = batch_data
                x = x.unsqueeze(2)  # Add channel dimension: [N, t, 1, 64, 64, 64]
                y = y.unsqueeze(1)  # Add channel dimension: [N, 1, 64, 64, 64]
            else:
                # Fallback for original format
                x, y = batch_data
            
            x, y = x.to(DEVICE, non_blocking=True), y.to(DEVICE, non_blocking=True)
            
            optimizer.zero_grad()
            pred = model(x)
            
            # Computing individual losses (same as original)
            l1 = F.l1_loss(y, pred) * L1_scale
            hist_l = hist_loss(y, pred, bins=32) * hist_scale
            
            mass_t = torch.pow(10, y).sum()
            mass_p = torch.pow(10, pred).sum()
            mass_error = torch.abs(mass_p - mass_t) / (mass_t + 1e-8)
            mass_l = torch.log1p(mass_error) * mass_scale
            
            spectral_l = compute_spectral_loss(y, pred) * spectral_scale
            
            with torch.no_grad():
                hd_true = torch.heaviside(y-torch.quantile(y.flatten(), 0.99), torch.tensor([1.0], device=DEVICE))
                hd_pred = torch.heaviside(pred-torch.quantile(pred.flatten(), 0.99), torch.tensor([1.0], device=DEVICE))
            hd_loss = F.l1_loss(hd_pred, hd_true)
            hd_l = hd_loss * hd_scale

            losses = torch.stack([l1, hist_l, mass_l, spectral_l, hd_l])
            if args.loss_type == 'static':
                total_loss = (static_w * losses).sum()
            else:
                w_t = torch.tensor(w, device=DEVICE)
                total_loss = (w_t * losses).sum()

            total_loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            
            # Accumulate losses
            batch_size = x.size(0)
            total_epoch_loss += total_loss.item() * batch_size
            sum_l1 += l1.item() * batch_size
            sum_hist += hist_l.item() * batch_size
            sum_mass += mass_l.item() * batch_size
            sum_spectral += spectral_l.item() * batch_size
            sum_hd += hd_l.item() * batch_size

            #timing printing
            batch_time_ema = print_progress(batch_idx, loader, epoch, start_epoch, args,
                                total_loss, batch_start, training_start,
                                batch_time_ema, ema_alpha)
            if batch_idx % 5 == 0:
                torch.cuda.empty_cache()

        N = len(dataset)
        avg_total = total_epoch_loss / N
        avg_l1 = sum_l1 / N
        avg_hist = sum_hist / N
        avg_mass = sum_mass / N
        avg_spectral = sum_spectral / N
        avg_hd = sum_hd / N
        
        if args.loss_type == 'dynamic':
            avg_losses = [avg_l1, avg_hist, avg_mass, avg_spectral, avg_hd]
            if epoch >= 2:
                r = [prev_losses[i]/prev2_losses[i] for i in range(n_losses)]
                T = args.DWA_temperature
                K = n_losses
                exp_r = [np.exp(r_i/T) for r_i in r]
                w = [(K * e) / sum(exp_r) for e in exp_r]
            else:
                w = [1.0]*n_losses
            prev2_losses = prev_losses
            prev_losses = avg_losses
            w_print = w
        else:
            w_print = static_w

        loss_history.append([float(avg_total), float(avg_l1), float(avg_hist), float(avg_mass), float(avg_spectral), float(avg_hd)])

        print(f"Epoch {epoch+1}: "
              f"Total={avg_total:.4f} | "
              f"L1={avg_l1:.4f} "
              f"HighDensity={avg_hd:.4f} "
              f"Mass={avg_mass:.4f}  Hist={avg_hist:.4f} "
              f"Spectral={avg_spectral:.4f} "
              f"Weights: {*w_print,}")
        print_memory_stats()
        # TensorBoard logging
        writer.add_scalar("Loss/Total", avg_total, epoch)
        writer.add_scalar("Loss/L1", avg_l1, epoch)
        writer.add_scalar("Loss/Hist", avg_hist, epoch)
        writer.add_scalar("Loss/Mass", avg_mass, epoch)
        writer.add_scalar("Loss/Spectral", avg_spectral, epoch)
        writer.add_scalar("Loss/HighDensity", avg_hd, epoch)

        if args.loss_type != 'static':
            curr_w = torch.softmax(loss_w, dim=0).detach().cpu().tolist()
            for i, name in enumerate(["L1", "Hist", "Mass", "Spectral", "HighDensity"]):
                writer.add_scalar(f"Weights/{name}", curr_w[i], epoch)

    writer.close()

    # Save model and plot
    torch.save({
        'epoch': args.epochs,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'loss': loss_history
    }, f"models/{args.run_name}_{MODELFILE}")
    
    if args.loss_type != 'static':
        torch.save(loss_w.detach().cpu(), f"models/plots/loss_data_{args.run_name}_loss.pt")
    
    import matplotlib.pyplot as plt
    loss_history = np.array(loss_history)
    plt.plot(loss_history[:,0], marker='o', label='Total Loss')
    plt.plot(loss_history[:,1], marker='.', label='L1 Loss')
    plt.plot(loss_history[:,2], marker='*', label='Histogram Loss')
    plt.plot(loss_history[:,3], marker='+', label='Mass Loss')
    plt.plot(loss_history[:,4], marker='x', label='Spectral Loss')
    plt.plot(loss_history[:,5], marker='1', label='High Density Loss')
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.xlabel("Epoch")
    plt.ylabel("Avg Weighted Loss")
    plt.title(f"{args.run_name} Training Loss")
    plt.savefig(f"models/plots/training_loss_{args.run_name}.png")
    plt.close()

    return min(loss_history[:,0])


def train_unet_with_dataset_optimized(dataset, args, argsGRU):
    '''
    Optimized training function with mixed precision and improved memory management.
    '''
    torch.cuda.empty_cache()
    torch.backends.cudnn.benchmark = False
    print(f"Training convGRUNet with {len(dataset)} samples, {args.run_name} run", flush=True)
    
    num_workers = args.num_workers if hasattr(args, 'num_workers') else min(8, os.cpu_count()//2)
    print(f"Using {num_workers} DataLoader workers", flush=True)
       
    loader = DataLoader(dataset,
                       batch_size=args.batch_size,
                       shuffle=True,
                       num_workers=num_workers,
                       pin_memory=True,
                       persistent_workers=(num_workers>0))
    
    # Initialize model, loss function, and optimizer
    model = hybridConvGRU3DNET(args, argsGRU).to(DEVICE, non_blocking=True)
    
    # For trainable weights
    n_losses = 5  # L1, hist, Mass, Spectral, High Density
    loss_w = torch.nn.Parameter(torch.ones(n_losses, device=DEVICE), requires_grad=True)
    optimizer = optim.AdamW(list(model.parameters()) + [loss_w], lr=args.lr, weight_decay=args.Adamw_weight_decay)
    writer = SummaryWriter(log_dir=f"models/logs/{args.run_name}")
    scaler = torch.amp.GradScaler(device=DEVICE)


    start_epoch = 0
    loss_history = []
    
    # Check if model already exists
    if os.path.exists(f"models/{args.run_name}_{MODELFILE}"):
        checkpoint = torch.load(f"models/{args.run_name}_{MODELFILE}", map_location=DEVICE, weights_only=False)
        model.load_state_dict(checkpoint['model_state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        start_epoch = checkpoint['epoch']
        loss_history = checkpoint['loss']
        print(f"Model loaded for training with {count_parameters(model):,} trainable parameters.", flush=True)
    else:
        print(f"Model initialized for training with {count_parameters(model):,} trainable parameters.", flush=True)
    
    # Loss scales
    L1_scale = 5.0
    hist_scale = 1.0
    mass_scale = 0.1
    spectral_scale = 0.5
    hd_scale = 1.0

    if args.loss_type == 'static':
        static_w = [1.0, 1.0, 1.0, 1.0, 1.0]
    if args.loss_type == 'dynamic':
        prev_losses = [1.0]*n_losses
        prev2_losses = [1.0]*n_losses
        w = [1.0]*n_losses

    # Training loop
    import time, datetime
    training_start = time.perf_counter()
    batch_time_ema = None
    ema_alpha = 0.15
    for epoch in range(start_epoch, args.epochs):
        model.train()
        total_epoch_loss = 0.0
        sum_l1, sum_hist, sum_mass, sum_spectral, sum_hd = 0.0, 0.0, 0.0, 0.0, 0.0
        torch.cuda.empty_cache()
        for batch_idx, batch_data in enumerate(loader):
            batch_start = time.perf_counter()
            if len(batch_data) == 3:
                x, y, labels = batch_data
                x = x.unsqueeze(2)  # Add channel dimension: [N, t, 1, 64, 64, 64]
                y = y.unsqueeze(1)  # Add channel dimension: [N, 1, 64, 64, 64]
            else:
                # Fallback for original format
                x, y = batch_data
            
            x, y = x.to(DEVICE, non_blocking=True), y.to(DEVICE, non_blocking=True)
            
            optimizer.zero_grad()
            with torch.amp.autocast(device_type=DEVICE.type, enabled=True):
                pred = model(x)
                # Computing individual losses (same as original)
                l1 = F.l1_loss(y, pred) * L1_scale
                hist_l = hist_loss(y, pred, bins=32) * hist_scale
                
                mass_t = torch.pow(10, y).sum()
                mass_p = torch.pow(10, pred).sum()
                mass_error = torch.abs(mass_p - mass_t) / (mass_t + 1e-8)
                mass_l = torch.log1p(mass_error) * mass_scale
                
                spectral_l = compute_spectral_loss(y, pred) * spectral_scale
                
                with torch.no_grad():
                    hd_true = torch.heaviside(y-torch.quantile(y.flatten(), 0.99), torch.tensor([1.0], device=DEVICE))
                    hd_pred = torch.heaviside(pred-torch.quantile(pred.flatten(), 0.99), torch.tensor([1.0], device=DEVICE))
                hd_loss = F.l1_loss(hd_pred, hd_true)
                hd_l = hd_loss * hd_scale

                losses = torch.stack([l1, hist_l, mass_l, spectral_l, hd_l])
                if args.loss_type == 'static':
                    total_loss = (static_w * losses).sum()
                else:
                    w_t = torch.tensor(w, device=DEVICE)
                    total_loss = (w_t * losses).sum()

            scaler.scale(total_loss).backward()
            scaler.unscale_(optimizer)  # if you clip grads
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            scaler.step(optimizer)
            scaler.update()
            
            # Accumulate losses
            batch_size = x.size(0)
            total_epoch_loss += total_loss.item() * batch_size
            sum_l1 += l1.item() * batch_size
            sum_hist += hist_l.item() * batch_size
            sum_mass += mass_l.item() * batch_size
            sum_spectral += spectral_l.item() * batch_size
            sum_hd += hd_l.item() * batch_size

            #timing printing
            batch_time_ema = print_progress(batch_idx, loader, epoch, start_epoch, args,
                                total_loss, batch_start, training_start,
                                batch_time_ema, ema_alpha)
            if batch_idx % 5 == 0:
                torch.cuda.empty_cache()

        N = len(dataset)
        avg_total = total_epoch_loss / N
        avg_l1 = sum_l1 / N
        avg_hist = sum_hist / N
        avg_mass = sum_mass / N
        avg_spectral = sum_spectral / N
        avg_hd = sum_hd / N
        
        if args.loss_type == 'dynamic':
            avg_losses = [avg_l1, avg_hist, avg_mass, avg_spectral, avg_hd]
            if epoch >= 2:
                r = [prev_losses[i]/prev2_losses[i] for i in range(n_losses)]
                T = args.DWA_temperature
                K = n_losses
                exp_r = [np.exp(r_i/T) for r_i in r]
                w = [(K * e) / sum(exp_r) for e in exp_r]
            else:
                w = [1.0]*n_losses
            prev2_losses = prev_losses
            prev_losses = avg_losses
            w_print = w
        else:
            w_print = static_w

        loss_history.append([float(avg_total), float(avg_l1), float(avg_hist), float(avg_mass), float(avg_spectral), float(avg_hd)])

        print(f"Epoch {epoch+1}: "
              f"Total={avg_total:.4f} | "
              f"L1={avg_l1:.4f} "
              f"HighDensity={avg_hd:.4f} "
              f"Mass={avg_mass:.4f}  Hist={avg_hist:.4f} "
              f"Spectral={avg_spectral:.4f} "
              f"Weights: {*w_print,}")
        print_memory_stats()
        # TensorBoard logging
        writer.add_scalar("Loss/Total", avg_total, epoch)
        writer.add_scalar("Loss/L1", avg_l1, epoch)
        writer.add_scalar("Loss/Hist", avg_hist, epoch)
        writer.add_scalar("Loss/Mass", avg_mass, epoch)
        writer.add_scalar("Loss/Spectral", avg_spectral, epoch)
        writer.add_scalar("Loss/HighDensity", avg_hd, epoch)

        if args.loss_type != 'static':
            curr_w = torch.softmax(loss_w, dim=0).detach().cpu().tolist()
            for i, name in enumerate(["L1", "Hist", "Mass", "Spectral", "HighDensity"]):
                writer.add_scalar(f"Weights/{name}", curr_w[i], epoch)

    writer.close()

    # Save model and plot
    torch.save({
        'epoch': args.epochs,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'loss': loss_history
    }, f"models/{args.run_name}_{MODELFILE}")
    
    if args.loss_type != 'static':
        torch.save(loss_w.detach().cpu(), f"models/plots/loss_data_{args.run_name}_loss.pt")
    
    import matplotlib.pyplot as plt
    loss_history = np.array(loss_history)
    plt.plot(loss_history[:,0], marker='o', label='Total Loss')
    plt.plot(loss_history[:,1], marker='.', label='L1 Loss')
    plt.plot(loss_history[:,2], marker='*', label='Histogram Loss')
    plt.plot(loss_history[:,3], marker='+', label='Mass Loss')
    plt.plot(loss_history[:,4], marker='x', label='Spectral Loss')
    plt.plot(loss_history[:,5], marker='1', label='High Density Loss')
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.xlabel("Epoch")
    plt.ylabel("Avg Weighted Loss")
    plt.title(f"{args.run_name} Training Loss")
    plt.savefig(f"models/plots/training_loss_{args.run_name}.png")
    plt.close()

    return min(loss_history[:,0])
