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
        x0_g = self.GRU3D1(x0)

        #First Encoder
        x1 = self.timeparser(x0, self.Conv2, doMaxPool = True)
        x1_g = self.GRU3D2(x1)

        if hasattr(self, 'Conv3'):
            #Second Encoder
            x2 = self.timeparser(x1, self.Conv3, doMaxPool = True)
            x2_g = self.GRU3D3(x2)
        if hasattr(self, 'Conv4'):
            #Third Encoder
            x3 = self.timeparser(x2, self.Conv4, doMaxPool = True)
            x3_g = self.GRU3D4(x3)
        if hasattr(self, 'Conv5'):
            #Fourth Encoder
            x4 = self.timeparser(x3, self.Conv5, doMaxPool = True)
            x4_g = self.GRU3D5(x4)

        all_oups = []
        for t in range(x.shape[1]):
            #decoding + concat path
            if hasattr(self, 'Conv5'):
                d3 = self.Up4(x4_g[:,t,:,:,:,:])
                a4 = self.Att4(g=d3,x=x3_g[:,t,:,:,:,:])
                d3 = torch.cat((a4,d3),dim=1)
                x3_g[:,t,:,:,:,:] = self.Up_conv4(d3)

            if hasattr(self, 'Conv4'):
                d2 = self.Up3(x3_g[:,t,:,:,:,:])
                a3 = self.Att3(g=d2,x=x2_g[:,t,:,:,:,:])
                d2 = torch.cat((a3,d2),dim=1)
                x2_g[:,t,:,:,:,:] = self.Up_conv3(d2)

            if hasattr(self, 'Conv3'):
                d1 = self.Up2(x2_g[:,t,:,:,:,:])
                a2 = self.Att2(g=d1,x=x1_g[:,t,:,:,:,:])
                d1 = torch.cat((a2,d1),dim=1)
                x1_g[:,t,:,:,:,:] = self.Up_conv2(d1)

            d0 = self.Up1(x1_g[:,t,:,:,:,:])
            a1 = self.Att1(g=d0,x=x0_g[:,t,:,:,:,:])
            d0 = torch.cat((a1,d0),dim=1)
            d0 = self.Up_conv1(d0)

            # Final 1x1 convolution to get the output channels
            all_oups.append(self.Conv_1x1(d0))
        
        out = torch.stack(all_oups, dim=1)
        return out


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
        out_arr = last_layer_output[:,-self.out_seq:,:,:,:,:] # shape (b,h,l)
        out_arr_i = []
        for out in range(self.out_seq):
            out_arr_i.append(self.out_cnn(out_arr[:,out,:,:,:,:])) # shape (b,t_out,c,l)
        out_arr = torch.stack(out_arr_i, dim=1)
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


