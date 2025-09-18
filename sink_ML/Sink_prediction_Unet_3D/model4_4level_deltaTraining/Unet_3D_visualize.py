# standard system modules
import os, sys
import h5py 
import argparse
# standard module for tabular data
import pandas as pd
import json

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
from torch.utils.data import DataLoader, TensorDataset
from torchvision.utils import save_image

from tqdm import tqdm

# set a seed to ensure reproducibility
seed = 128
rnd  = np.random.RandomState(seed)


DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
#DEVICE = torch.device('cpu')

print(f'Available device: {str(DEVICE):4s}', flush=True)

def count_parameters(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


import torch.nn as nn

class AttentionGate(nn.Module):
    def __init__(self, F_g, F_l, F_int):
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
        self.relu = nn.ReLU(inplace=True)

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
    def __init__(self,ch_in,ch_out):
        super().__init__()
        self.conv = nn.Sequential(
            nn.Conv3d(ch_in, ch_out, kernel_size=3,stride=1,padding=1,bias=True, padding_mode='circular'),
            nn.BatchNorm3d(ch_out),
            nn.ReLU(inplace=True),
            nn.Conv3d(ch_out, ch_out, kernel_size=3,stride=1,padding=1,bias=True, padding_mode='circular'),
            nn.BatchNorm3d(ch_out),
            nn.ReLU(inplace=True)
        )

    def forward(self,x):
        x = self.conv(x)
        return x

class Upsampling_Block(nn.Module):
    def __init__(self,ch_in,ch_out):
        super().__init__()
        self.up = nn.Sequential(
            nn.Upsample(scale_factor=2),
            nn.Conv3d(ch_in,ch_out,kernel_size=3,stride=1,padding=1,bias=True, padding_mode='circular'),
		    nn.BatchNorm3d(ch_out),
			nn.ReLU(inplace=True)
        )

    def forward(self,x):
        x = self.up(x)
        return x


class UNet3D(nn.Module):
    def __init__(self, in_ch, out_ch, base_f=32):
        """
        A 3D U-Net with attention gates, 4 level deep
        Args:
            in_ch: number of input channels (e.g., 1 for grayscale)
            out_ch: number of output channels (e.g., 1 for grayscale)
            base_f: base number of filters in the first layer
        """
        super().__init__()

        self.Maxpool = nn.MaxPool3d(kernel_size=2,stride=2)

        self.Conv1 = Convolution_Block(ch_in=in_ch,ch_out=base_f)
        self.Conv2 = Convolution_Block(ch_in=base_f,ch_out=base_f*2)
        self.Conv3 = Convolution_Block(ch_in=base_f*2,ch_out=base_f*4)
        self.Conv4 = Convolution_Block(ch_in=base_f*4,ch_out=base_f*8)
        '''
        To add a deeper level with attention, add levels like this lol
        self.Conv5 = Convolution_Block(ch_in=base_f*8,ch_out=base_f*16)

        self.Up4 = Upsampling_Block(ch_in=base_f*16,ch_out=base_f*8)
        self.Att4 = AttentionGate(F_g=base_f*8,F_l=base_f*8,F_int=base_f*4)
        self.Up_conv4 = Convolution_Block(ch_in=base_f*16, ch_out=base_f*8)
        '''

        self.Up3 = Upsampling_Block(ch_in=base_f*8,ch_out=base_f*4)
        self.Att3 = AttentionGate(F_g=base_f*4,F_l=base_f*4,F_int=base_f*2)
        self.Up_conv3 = Convolution_Block(ch_in=base_f*8, ch_out=base_f*4)

        self.Up2 = Upsampling_Block(ch_in=base_f*4,ch_out=base_f*2)
        self.Att2 = AttentionGate(F_g=base_f*2,F_l=base_f*2,F_int=base_f)
        self.Up_conv2 = Convolution_Block(ch_in=base_f*4, ch_out=base_f*2)

        self.Up1 = Upsampling_Block(ch_in=base_f*2,ch_out=base_f)
        self.Att1 = AttentionGate(F_g=base_f,F_l=base_f,F_int=base_f//2)
        self.Up_conv1 = Convolution_Block(ch_in=base_f*2, ch_out=base_f)

        self.Conv_1x1 = nn.Conv3d(base_f, out_ch, kernel_size=1, stride=1, padding=0)

    def forward(self,x):
        # encoding path
        x0 = self.Conv1(x)

        #First Encoder
        x1 = self.Maxpool(x0)
        x1 = self.Conv2(x1)

        #Second Encoder
        x2 = self.Maxpool(x1)
        x2 = self.Conv3(x2)

        #Third Encoder
        x3 = self.Maxpool(x2)
        x3 = self.Conv4(x3)


        #decoding + concat path
        d2 = self.Up3(x3)
        a3 = self.Att3(g=d2,x=x2)
        d2 = torch.cat((a3,d2),dim=1)
        d2 = self.Up_conv3(d2)

        d1 = self.Up2(d2)
        a2 = self.Att2(g=d1,x=x1)
        d1 = torch.cat((a2,d1),dim=1)
        d1 = self.Up_conv2(d1)

        d0 = self.Up1(d1)
        a1 = self.Att1(g=d0,x=x0)
        d0 = torch.cat((a1,d0),dim=1)
        d0 = self.Up_conv1(d0)

        # Final 1x1 convolution to get the output channels
        return self.Conv_1x1(d0)

from torchview import draw_graph

parser = argparse.ArgumentParser()
args = parser.parse_args(args=[])
args.batch_size = 8
args.image_size = IMAGESIZE
args.device = DEVICE
args.run_name = "Unet_3D_HD_Delta_scaled"
args.epochs = 30
args.lr = 3e-4
args.loss_type = 'dynamic' # 'static' or 'dynamic'


model = UNet3D(in_ch=1, out_ch=1)
batch_size = 8
# device='meta' -> no memory is consumed for visualization
model_graph = draw_graph(model, input_size=(batch_size,1,64,64,64), device='meta', expand_nested=True, depth=1)
import graphviz
graph_data = str(model_graph.visual_graph)

fie_ext = 'pdf'
my_graph = graphviz.Source(graph_data)
my_graph.render(args.run_name,format=fie_ext, view=False)