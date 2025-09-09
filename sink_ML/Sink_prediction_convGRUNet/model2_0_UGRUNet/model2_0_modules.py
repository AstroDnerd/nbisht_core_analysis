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
    def __init__(self, in_ch, out_ch, args, argsGRU):
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
        self.GRU3D1 = ConvGRUNet3D(input_channels=args.base_f, output_channels =args.base_f, image_size=args.image_size, in_seq = argsGRU.in_seq,
                                      out_seq = argsGRU.out_seq, num_layers = argsGRU.num_layers, kernels = argsGRU.kernels, hidden_dims = argsGRU.hidden_dims, 
                                      bias = argsGRU.bias)
        self.Conv2 = Convolution_Block(ch_in=args.base_f,ch_out=args.base_f*2, args=args)
        self.GRU3D2 = ConvGRUNet3D(input_channels=args.base_f*2, output_channels =args.base_f*2, image_size=args.image_size//2, in_seq = argsGRU.in_seq,
                                      out_seq = argsGRU.out_seq, num_layers = argsGRU.num_layers, kernels = argsGRU.kernels, hidden_dims = argsGRU.hidden_dims, 
                                      bias = argsGRU.bias)
        if args.depth >= 3:
            self.Conv3 = Convolution_Block(ch_in=args.base_f*2,ch_out=args.base_f*4, args=args)
            self.GRU3D3 = ConvGRUNet3D(input_channels=args.base_f*4, output_channels =args.base_f*4, image_size=args.image_size//4, in_seq = argsGRU.in_seq,
                                        out_seq = argsGRU.out_seq, num_layers = argsGRU.num_layers, kernels = argsGRU.kernels, hidden_dims = argsGRU.hidden_dims, 
                                        bias = argsGRU.bias)
        if args.depth >= 4:
            self.Conv4 = Convolution_Block(ch_in=args.base_f*4,ch_out=args.base_f*8, args=args)
            self.GRU3D4 = ConvGRUNet3D(input_channels=args.base_f*8, output_channels =args.base_f*8, image_size=args.image_size//8, in_seq = argsGRU.in_seq,
                                        out_seq = argsGRU.out_seq, num_layers = argsGRU.num_layers, kernels = argsGRU.kernels, hidden_dims = argsGRU.hidden_dims, 
                                        bias = argsGRU.bias)
        if args.depth >= 5:
            self.Conv5 = Convolution_Block(ch_in=args.base_f*8,ch_out=args.base_f*16, args=args)
            self.GRU3D5 = ConvGRUNet3D(input_channels=args.base_f*16, output_channels =args.base_f*16, image_size=args.image_size//16, in_seq = argsGRU.in_seq,
                                        out_seq = argsGRU.out_seq, num_layers = argsGRU.num_layers, kernels = argsGRU.kernels, hidden_dims = argsGRU.hidden_dims, 
                                        bias = argsGRU.bias)

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

    def timeparser(self, xinput, unetfunc, doMaxPool = True):
        #xinput of form (b,t,c,l,h,w)
        all_inps = []
        for t in range(xinput.shape[1]):
            xi = xinput[:,t,:,:,:,:] #xi is (b,c,l,h,w)
            if doMaxPool:
                xi = self.Maxpool(xi)
            xj = unetfunc(xi)
            all_inps.append(xj)
        xj_stacked = torch.stack(all_inps, dim=1)
        return xj_stacked

    def forward(self,x):
        # encoding path

        x0 = self.timeparser(x, self.Conv1, doMaxPool = False)
        x0_g = self.GRU3D1(x0)[:,-1,:,:,:,:]

        #First Encoder
        x1 = self.timeparser(x0, self.Conv2, doMaxPool = True)
        x1_g = self.GRU3D2(x1)[:,-1,:,:,:,:]

        if hasattr(self, 'Conv3'):
            #Second Encoder
            x2 = self.timeparser(x1, self.Conv3, doMaxPool = True)
            x2_g = self.GRU3D3(x2)[:,-1,:,:,:,:]
        if hasattr(self, 'Conv4'):
            #Third Encoder
            x3 = self.timeparser(x2, self.Conv4, doMaxPool = True)
            x3_g = self.GRU3D4(x3)[:,-1,:,:,:,:]
        if hasattr(self, 'Conv5'):
            #Fourth Encoder
            x4 = self.timeparser(x3, self.Conv5, doMaxPool = True)
            x4_g = self.GRU3D5(x4)[:,-1,:,:,:,:]


        #decoding + concat path
        if hasattr(self, 'Conv5'):
            d3 = self.Up4(x4_g)
            a4 = self.Att4(g=d3,x=x3_g)
            d3 = torch.cat((a4,d3),dim=1)
            x3_g = self.Up_conv4(d3)

        if hasattr(self, 'Conv4'):
            d2 = self.Up3(x3_g)
            a3 = self.Att3(g=d2,x=x2_g)
            d2 = torch.cat((a3,d2),dim=1)
            x2_g = self.Up_conv3(d2)

        if hasattr(self, 'Conv3'):
            d1 = self.Up2(x2_g)
            a2 = self.Att2(g=d1,x=x1_g)
            d1 = torch.cat((a2,d1),dim=1)
            x1_g = self.Up_conv2(d1)

        d0 = self.Up1(x1_g)
        a1 = self.Att1(g=d0,x=x0_g)
        d0 = torch.cat((a1,d0),dim=1)
        d0 = self.Up_conv1(d0)

        # Final 1x1 convolution to get the output channels
        return self.Conv_1x1(d0)


class hybridUGRUNET(nn.Module):
    def __init__(self, args, argsGRU):
        """
        Initialize the hybridUGRUNET model
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

        self.UNet3D = UNet3D(in_ch=args.in_channel, out_ch=args.out_channel, args=args, argsGRU = argsGRU)

    def forward(self, input_tensor):
        """
        :param input_tensor: (b, t_in, c, l)
        :param hidden_state:
        :return: output_tensor: (b, t_out, c, l)
        """
        output = self.UNet3D(input_tensor)  #(b, t_out, c, l, h ,w)
        return output
        

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

def return_total_loss(y, pred, L1_scale = 5.0, hist_scale = 1.0, mass_scale = 0.1, spectral_scale = 0.5, hd_scale = 1.0):
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

    return l1, hist_l, mass_l, spectral_l, hd_l

def check_training_phase(epoch):
    #EPOCH 0 losses for 16^3 sim
    #L1=1.4236 HighDensity=0.0164 Mass=0.0138  Hist=0.0064 Spectral=0.3902
    L1_scale, hist_scale, mass_scale, spectral_scale, hd_scale = 5.0, 1.0, 0.1, 0.5, 1.0
    epoch_ranges = [50,70,80,90,100]
    if 0<=epoch<=epoch_ranges[0]:
        L1_scale = 10*L1_scale
    elif epoch_ranges[0]<epoch<=epoch_ranges[1]:
        mass_scale = 100*mass_scale
    elif epoch_ranges[1]<epoch<=epoch_ranges[2]:
        spectral_scale = 5*spectral_scale
    elif epoch_ranges[2]<epoch<=epoch_ranges[3]:
        hist_scale = 100*hist_scale
    elif epoch_ranges[3]<epoch<=epoch_ranges[4]:
        hd_scale = 100*hd_scale
    return L1_scale, hist_scale, mass_scale,spectral_scale,hd_scale


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
    model = hybridUGRUNET(args, argsGRU).to(DEVICE, non_blocking=True)
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

def train_unet(input_arr, output_arr, labels, args, argsGRU, validation_input = None, validation_output = None):
    print(f"Training convGRUNet with {len(input_arr)} samples, {args.run_name} run", flush=True)
    X = torch.stack(input_arr, dim=0).float().unsqueeze(2)  # [N, t, 1, 64, 64, 64]
    Y = torch.stack(output_arr, dim=0).float().unsqueeze(1) # [N, 1, 64, 64, 64]

    if validation_input!=None:
        #validation
        val_x_arr = torch.stack(validation_input, dim=0).float().unsqueeze(2)  # [N, t, 1, 16, 16, 16]
        val_y_arr = torch.stack(validation_output, dim=0).float().unsqueeze(1)  # [N, t, 1, 16, 16, 16]

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
    model = hybridUGRUNET(args, argsGRU).to(DEVICE, non_blocking=True)
    # For trainable weights
    n_losses = 5  # L1, hist, Mass, Spectral, High Density
    loss_w = torch.nn.Parameter(torch.ones(n_losses, device=DEVICE), requires_grad=True)
    optimizer = optim.AdamW(list(model.parameters()) + [loss_w], lr=args.lr, weight_decay=args.Adamw_weight_decay)
    LRSched = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=2, threshold=0.0001, threshold_mode='rel', cooldown=0, min_lr=0, eps=1e-08)
    writer = SummaryWriter(log_dir=f"models/logs/{args.run_name}")

    start_epoch = 0
    loss_history = []
    val_loss_history = []
    #check if model already exists
    if os.path.exists(f"models/{args.run_name}_{MODELFILE}"):
        checkpoint = torch.load(f"models/{args.run_name}_{MODELFILE}", map_location=DEVICE, weights_only=False)
        model.load_state_dict(checkpoint['model_state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        LRSched.load_state_dict(checkpoint['scheduler_state_dict'])
        start_epoch = checkpoint['epoch']
        loss_history = checkpoint['loss']
        val_loss_history = checkpoint['val_loss']
        print(f"Model loaded for training with {count_parameters(model):,} trainable parameters.", flush=True)
    else:
        print(f"Model initialized for training with {count_parameters(model):,} trainable parameters.", flush=True)
    
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
            l1, hist_l, mass_l, spectral_l, hd_l = return_total_loss(y, pred)

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
        
        last_lr = LRSched.get_last_lr()
        LRSched.step(avg_total)

        loss_history.append([float(avg_total), float(avg_l1), float(avg_hist), float(avg_mass), float(avg_spectral), float(avg_hd)])

        # log to console
        print(f"Epoch {epoch+1}: "
              f"Total={avg_total:.4f} | "
              f"L1={avg_l1:.4f} "
              f"HighDensity={avg_hd:.4f} "
              f"Mass={avg_mass:.4f}  Hist={avg_hist:.4f} "
              f"Spectral={avg_spectral:.4f} "
              f"LR={*last_lr,} "
              f"Weights: {*w_print,}", flush=True)

        # log to TensorBoard
        writer.add_scalar("Loss/Total",    avg_total, epoch)
        writer.add_scalar("Loss/L1",       avg_l1,    epoch)
        writer.add_scalar("Loss/Hist",     avg_hist,   epoch)
        writer.add_scalar("Loss/Mass",     avg_mass,  epoch)
        writer.add_scalar("Loss/Spectral", avg_spectral, epoch)
        writer.add_scalar("Loss/HighDensity", avg_hd, epoch)
        writer.add_scalar("LearningRate", last_lr[-1], epoch)

        if args.loss_type!='static':
            # log the learned weights (after softmax)
            curr_w = torch.softmax(loss_w, dim=0).detach().cpu().tolist()
            for i, name in enumerate(["L1","Hist","Mass", "Spectral", "HighDensity"]):
                writer.add_scalar(f"Weights/{name}", curr_w[i], epoch)
        
        if epoch%5==0 and validation_input!=None:
            val_total, val_l1, val_hist, val_mass, val_spectral, val_hd = 0, 0, 0, 0, 0, 0
            for i in range(len(validation_input)):
                val_x = val_x_arr[i].unsqueeze(0)
                val_y = val_y_arr[i].unsqueeze(0)
                with torch.no_grad():
                    pred = model(val_x)
                    l1, hist_l, mass_l, spectral_l, hd_l = return_total_loss(val_y, pred)
                    losses = torch.stack([l1, hist_l, mass_l, spectral_l, hd_l])
                    w_t = torch.tensor(w_print)
                    total_loss = (w_t * losses).sum()
                    val_total += total_loss.item()
                    val_l1     += l1.item()
                    val_hist    += hist_l.item()
                    val_mass   += mass_l.item()
                    val_spectral   += spectral_l.item()
                    val_hd   += hd_l.item()
            N = len(validation_input)
            avg_total = val_total / N
            avg_l1    = val_l1     / N
            avg_hist  = val_hist   / N
            avg_mass  = val_mass   / N
            avg_spectral = val_spectral / N
            avg_hd = val_hd / N

            # log to console
            print(f"Validation! Epoch {epoch+1}: "
                f"Total={avg_total:.4f} | "
                f"L1={avg_l1:.4f} "
                f"HighDensity={avg_hd:.4f} "
                f"Mass={avg_mass:.4f}  Hist={avg_hist:.4f} "
                f"Spectral={avg_spectral:.4f} "
                f"LR={*last_lr,} "
                f"Weights: {*w_print,}", flush=True)
            val_loss_history.append([float(avg_total), float(avg_l1), float(avg_hist), float(avg_mass), float(avg_spectral), float(avg_hd)])

    writer.close()

    # === Save & plot ===
    torch.save({'epoch': args.epochs, 'model_state_dict': model.state_dict(), 'optimizer_state_dict': optimizer.state_dict(), 
                'scheduler_state_dict': LRSched.state_dict(), 'loss': loss_history, 'val_loss': val_loss_history}, f"models/{args.run_name}_{MODELFILE}")
    if args.loss_type != 'static':
        torch.save(loss_w.detach().cpu(), f"models/plots/loss_data_{args.run_name}_loss.pt")
    import matplotlib.pyplot as plt
    loss_history = np.array(loss_history)
    val_loss_history = np.array(val_loss_history)
    plt.figure()
    plt.plot(loss_history[:,0], marker='o', label='Total Loss')
    plt.plot(loss_history[:,1], marker='.', label='L1 Loss')
    plt.plot(loss_history[:,2], marker='*', label='Histogram Loss')
    plt.plot(loss_history[:,3], marker='+', label='Mass Loss')
    plt.plot(loss_history[:,4], marker='x', label='Spectral Loss')
    plt.plot(loss_history[:,5], marker='1', label='High Density Loss')
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

    plt.figure()
    plt.plot(val_loss_history[:,0], marker='o', label='Total Loss')
    plt.plot(val_loss_history[:,1], marker='.', label='L1 Loss')
    plt.plot(val_loss_history[:,2], marker='*', label='Histogram Loss')
    plt.plot(val_loss_history[:,3], marker='+', label='Mass Loss')
    plt.plot(val_loss_history[:,4], marker='x', label='Spectral Loss')
    plt.plot(val_loss_history[:,5], marker='1', label='High Density Loss')
    plt.plot(val_loss_history[:,5], marker='1', label='High Density Loss')
    plt.yscale('log')  # log scale for better visibility of loss changes
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.xlabel("Every 5th Epoch")
    plt.ylabel("Avg Weighted Loss")
    plt.title(f"{args.run_name} Validation Training Loss")
    plt.savefig(f"models/plots/validation_training_loss_{args.run_name}.png")
    plt.close()

    return min(loss_history[:,0])  # return the minimum loss for this run


def train_unet_training_phase(input_arr, output_arr, labels, args, argsGRU, validation_input = None, validation_output = None):
    print(f"Training convGRUNet with {len(input_arr)} samples, {args.run_name} run", flush=True)
    X = torch.stack(input_arr, dim=0).float().unsqueeze(2)  # [N, t, 1, 64, 64, 64]
    Y = torch.stack(output_arr, dim=0).float().unsqueeze(1) # [N, 1, 64, 64, 64]

    if validation_input!=None:
        #validation
        val_x_arr = torch.stack(validation_input, dim=0).float().unsqueeze(2)  # [N, t, 1, 64, 64, 64]
        val_y_arr = torch.stack(validation_output, dim=0).float().unsqueeze(1)  # [N, t, 1, 64, 64, 64]

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
    model = hybridUGRUNET(args, argsGRU).to(DEVICE, non_blocking=True)
    # For trainable weights
    n_losses = 5  # L1, hist, Mass, Spectral, High Density
    loss_w = torch.nn.Parameter(torch.ones(n_losses, device=DEVICE), requires_grad=True)
    optimizer = optim.AdamW(list(model.parameters()) + [loss_w], lr=args.lr, weight_decay=args.Adamw_weight_decay)
    LRSched = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=2, threshold=0.0001, threshold_mode='rel', cooldown=0, min_lr=0, eps=1e-08)
    writer = SummaryWriter(log_dir=f"models/logs/{args.run_name}")

    start_epoch = 0
    loss_history = []
    val_loss_history = []
    #check if model already exists
    if os.path.exists(f"models/{args.run_name}_{MODELFILE}"):
        checkpoint = torch.load(f"models/{args.run_name}_{MODELFILE}", map_location=DEVICE, weights_only=False)
        model.load_state_dict(checkpoint['model_state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        LRSched.load_state_dict(checkpoint['scheduler_state_dict'])
        start_epoch = checkpoint['epoch']
        loss_history = checkpoint['loss']
        val_loss_history = checkpoint['val_loss']
        print(f"Model loaded for training with {count_parameters(model):,} trainable parameters.", flush=True)
    else:
        print(f"Model initialized for training with {count_parameters(model):,} trainable parameters.", flush=True)
    
    if args.loss_type == 'static':
        static_w = [1.0, 1.0, 1.0, 1.0, 1.0]  # initial static weights for L1, hist, Mass, Spectral, High Density
    if args.loss_type == 'dynamic':
        prev_losses = [1.0]*n_losses    # dummy for epoch 0
        prev2_losses = [1.0]*n_losses   # dummy for epoch -1
        w = [1.0]*n_losses

    # === Training loop ===
    print("TRAINING PHASE 1!!!", flush=True)
    tp=2
    for epoch in range(start_epoch, args.epochs):
        model.train()
        total_epoch_loss = 0.0
        sum_l1, sum_hist, sum_mass, sum_spectral, sum_hd = 0.0, 0.0, 0.0, 0.0, 0.0
        #Training phases check epoch condition
        l1_s, hist_s, mass_s, spec_s, hd_s = check_training_phase(epoch)


        for x, y in loader:
            x, y = x.to(DEVICE, non_blocking=True), y.to(DEVICE, non_blocking=True)
            optimizer.zero_grad()
            pred = model(x)
            l1, hist_l, mass_l, spectral_l, hd_l = return_total_loss(y, pred, L1_scale = l1_s, hist_scale = hist_s, mass_scale = mass_s,
                                                                      spectral_scale = spec_s, hd_scale = hd_s)

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
        
        last_lr = LRSched.get_last_lr()
        LRSched.step(avg_total)

        loss_history.append([float(avg_total), float(avg_l1), float(avg_hist), float(avg_mass), float(avg_spectral), float(avg_hd)])

        # log to console
        print(f"Epoch {epoch+1}: "
              f"Total={avg_total:.4f} | "
              f"L1={avg_l1:.4f} "
              f"HighDensity={avg_hd:.4f} "
              f"Mass={avg_mass:.4f}  Hist={avg_hist:.4f} "
              f"Spectral={avg_spectral:.4f} "
              f"LR={*last_lr,} "
              f"Weights: {*w_print,}", flush=True)

        # log to TensorBoard
        writer.add_scalar("Loss/Total",    avg_total, epoch)
        writer.add_scalar("Loss/L1",       avg_l1,    epoch)
        writer.add_scalar("Loss/Hist",     avg_hist,   epoch)
        writer.add_scalar("Loss/Mass",     avg_mass,  epoch)
        writer.add_scalar("Loss/Spectral", avg_spectral, epoch)
        writer.add_scalar("Loss/HighDensity", avg_hd, epoch)
        writer.add_scalar("LearningRate", last_lr[-1], epoch)

        if args.loss_type!='static':
            # log the learned weights (after softmax)
            curr_w = torch.softmax(loss_w, dim=0).detach().cpu().tolist()
            for i, name in enumerate(["L1","Hist","Mass", "Spectral", "HighDensity"]):
                writer.add_scalar(f"Weights/{name}", curr_w[i], epoch)
        
        if epoch%5==0 and validation_input!=None:
            val_total, val_l1, val_hist, val_mass, val_spectral, val_hd = 0, 0, 0, 0, 0, 0
            for i in range(len(validation_input)):
                val_x = val_x_arr[i].unsqueeze(0)
                val_y = val_y_arr[i].unsqueeze(0)
                with torch.no_grad():
                    pred = model(val_x)
                    l1, hist_l, mass_l, spectral_l, hd_l = return_total_loss(val_y, pred, L1_scale = l1_s, hist_scale = hist_s, mass_scale = mass_s,
                                                                      spectral_scale = spec_s, hd_scale = hd_s)
                    losses = torch.stack([l1, hist_l, mass_l, spectral_l, hd_l])
                    w_t = torch.tensor(w_print)
                    total_loss = (w_t * losses).sum()
                    val_total += total_loss.item()
                    val_l1     += l1.item()
                    val_hist    += hist_l.item()
                    val_mass   += mass_l.item()
                    val_spectral   += spectral_l.item()
                    val_hd   += hd_l.item()
            N = len(validation_input)
            avg_total = val_total / N
            avg_l1    = val_l1     / N
            avg_hist  = val_hist   / N
            avg_mass  = val_mass   / N
            avg_spectral = val_spectral / N
            avg_hd = val_hd / N

            # log to console
            print(f"Validation! Epoch {epoch+1}: "
                f"Total={avg_total:.4f} | "
                f"L1={avg_l1:.4f} "
                f"HighDensity={avg_hd:.4f} "
                f"Mass={avg_mass:.4f}  Hist={avg_hist:.4f} "
                f"Spectral={avg_spectral:.4f} "
                f"LR={*last_lr,} "
                f"Weights: {*w_print,}", flush=True)
            val_loss_history.append([float(avg_total), float(avg_l1), float(avg_hist), float(avg_mass), float(avg_spectral), float(avg_hd)])
        
        epoch_ranges = [30,50,70,80,90]
        if epoch in epoch_ranges:
            if epoch!=epoch_ranges[0]:
                print("TRAINING PHASE "+str(tp)+"!!!", flush=True)
                tp+=1
                for g in optimizer.param_groups:
                    g['lr'] = args.lr
            torch.save({'epoch': args.epochs, 'model_state_dict': model.state_dict(), 'optimizer_state_dict': optimizer.state_dict(), 
                'scheduler_state_dict': LRSched.state_dict(), 'loss': loss_history, 'val_loss': val_loss_history}, f"models/{args.run_name}_{MODELFILE}")
            torch.save(loss_w.detach().cpu(), f"models/plots/loss_data_{args.run_name}_loss.pt")


    writer.close()

    # === Save & plot ===
    torch.save({'epoch': args.epochs, 'model_state_dict': model.state_dict(), 'optimizer_state_dict': optimizer.state_dict(), 
                'scheduler_state_dict': LRSched.state_dict(), 'loss': loss_history, 'val_loss': val_loss_history}, f"models/{args.run_name}_{MODELFILE}")
    if args.loss_type != 'static':
        torch.save(loss_w.detach().cpu(), f"models/plots/loss_data_{args.run_name}_loss.pt")
    import matplotlib.pyplot as plt
    loss_history = np.array(loss_history)
    val_loss_history = np.array(val_loss_history)
    plt.figure()
    plt.plot(loss_history[:,0], marker='o', label='Total Loss')
    plt.plot(loss_history[:,1], marker='.', label='L1 Loss')
    plt.plot(loss_history[:,2], marker='*', label='Histogram Loss')
    plt.plot(loss_history[:,3], marker='+', label='Mass Loss')
    plt.plot(loss_history[:,4], marker='x', label='Spectral Loss')
    plt.plot(loss_history[:,5], marker='1', label='High Density Loss')
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

    plt.figure()
    plt.plot(val_loss_history[:,0], marker='o', label='Total Loss')
    plt.plot(val_loss_history[:,1], marker='.', label='L1 Loss')
    plt.plot(val_loss_history[:,2], marker='*', label='Histogram Loss')
    plt.plot(val_loss_history[:,3], marker='+', label='Mass Loss')
    plt.plot(val_loss_history[:,4], marker='x', label='Spectral Loss')
    plt.plot(val_loss_history[:,5], marker='1', label='High Density Loss')
    plt.plot(val_loss_history[:,5], marker='1', label='High Density Loss')
    plt.yscale('log')  # log scale for better visibility of loss changes
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.xlabel("Every 5th Epoch")
    plt.ylabel("Avg Weighted Loss")
    plt.title(f"{args.run_name} Validation Training Loss")
    plt.savefig(f"models/plots/validation_training_loss_{args.run_name}.png")
    plt.close()

    return min(loss_history[:,0])  # return the minimum loss for this run
