import os
import copy
import numpy as np
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
import torch.nn.functional as F
from tqdm import tqdm
from torch import optim
import torchvision
import json
from torchvision.utils import save_image
import logging
from torch.utils.tensorboard import SummaryWriter

logging.basicConfig(format="%(asctime)s - %(levelname)s: %(message)s", level=logging.INFO, datefmt="%I:%M:%S")

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

def get_images(input_dir):
    output_dic = {}
    for filename in os.listdir(input_dir):
        infile = open(os.path.join(input_dir, filename), 'r')
        i_file = json.load(infile)
        num = filename.split('_')[1]
        img_dic = i_file['img']
        label = i_file['label']
        output_dic[num] = {'img': img_dic, 'label': label}
        infile.close()
    return output_dic


class EMA:
    def __init__(self, beta):
        super().__init__()
        self.beta = beta
        self.step = 0

    def update_model_average(self, ma_model, current_model):
        for current_params, ma_params in zip(current_model.parameters(), ma_model.parameters()):
            old_weight, up_weight = ma_params.data, current_params.data
            ma_params.data = self.update_average(old_weight, up_weight)

    def update_average(self, old, new):
        if old is None:
            return new
        return old * self.beta + (1 - self.beta) * new

    def step_ema(self, ema_model, model, step_start_ema=2000):
        if self.step < step_start_ema:
            self.reset_parameters(ema_model, model)
            self.step += 1
            return
        self.update_model_average(ema_model, model)
        self.step += 1

    def reset_parameters(self, ema_model, model):
        ema_model.load_state_dict(model.state_dict())



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

