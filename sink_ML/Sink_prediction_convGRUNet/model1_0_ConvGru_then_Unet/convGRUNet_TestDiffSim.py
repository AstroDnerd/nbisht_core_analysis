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

IMAGEINPUT = '/home/x-nbisht1/scratch/projects/128/B2/datasets/Unet/u502_Images_Input_3D/'
IMAGEOUTPUT = '/home/x-nbisht1/scratch/projects/128/B2/datasets/Unet/u502_Images_Output_3D/'
MODELFILE = 'nnmodel.dict'

IMAGESIZE = 64
FRAMES = np.arange(0,90, 1)
FRAME_DIFF = 30

DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
#DEVICE = torch.device('cpu')

print(f'Available device: {str(DEVICE):4s}', flush=True)

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
    for filename in os.listdir(input_dir)[::10]:  # Load every 100 file to reduce memory usage
        infile = open(os.path.join(input_dir, filename), 'r')
        i_file = json.load(infile)
        num = filename.split('_')[1]
        img_dic = i_file['img']
        label = i_file['label']
        output_dic[num] = {'img': img_dic, 'label': label}
        infile.close()
    return output_dic


print('Loading images!', flush=True)
input_image_dic = get_images(IMAGEINPUT)
output_image_dic = get_images(IMAGEOUTPUT)
input_image_arr = []
output_image_arr = []
label_arr = []
for key in input_image_dic.keys():
    input_image_arr.append(torch.tensor(input_image_dic[key]['img']))
    output_image_arr.append(torch.tensor(output_image_dic[key]['img']))
    label_arr.append(input_image_dic[key]['label'])


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

def plot_comparison(input_arr, pred_output_arr, act_output_arr, label, loss_dict, save_filename='xy_data.png', dont_show = False):
    input_str = 'INPUT IMAGE:\n'
    pred_str = 'PREDICTED IMAGE:\n'
    act_str = 'ACTUAL IMAGE:\n'
    # Calculate metrics

    input_L1 = F.l1_loss(torch.from_numpy(act_output_arr), torch.from_numpy(input_arr))
    loss_dict['l1']['input'].append(input_L1)
    input_str+= 'L1: '+'{:0.5f}'.format(input_L1)+'\n'
    pred_L1 = F.l1_loss(torch.from_numpy(act_output_arr), torch.from_numpy(pred_output_arr))
    loss_dict['l1']['pred'].append(pred_L1)
    pred_str+= 'L1: '+'{:0.5f}'.format(pred_L1)+'\n'
    actual_L1 = F.l1_loss(torch.from_numpy(act_output_arr), torch.from_numpy(act_output_arr))
    act_str+= 'L1: '+'{:0.5f}'.format(actual_L1)+'\n'

    input_hist_loss = hist_loss(torch.from_numpy(act_output_arr), torch.from_numpy(input_arr), bins=32)
    loss_dict['hist']['input'].append(input_hist_loss)
    input_str+= 'HistLoss: '+'{:0.5f}'.format(input_hist_loss)+'\n'
    pred_hist_loss = hist_loss(torch.from_numpy(act_output_arr), torch.from_numpy(pred_output_arr), bins=32)
    loss_dict['hist']['pred'].append(pred_hist_loss)
    pred_str+= 'HistLoss: '+'{:0.5f}'.format(pred_hist_loss)+'\n'
    actual_hist_loss = hist_loss(torch.from_numpy(act_output_arr), torch.from_numpy(act_output_arr), bins=32)
    act_str+= 'HistLoss: '+'{:0.5f}'.format(actual_hist_loss)+'\n'

    inversed_input_arr = img_inverse_transform(input_arr)
    inversed_pred_output_arr = img_inverse_transform(pred_output_arr)
    inversed_act_output_arr = img_inverse_transform(act_output_arr)
    input_total_density = np.sum(inversed_input_arr)
    pred_total_density = np.sum(inversed_pred_output_arr)
    act_total_density = np.sum(inversed_act_output_arr)
    input_mass_loss = abs(act_total_density - input_total_density) / (act_total_density+ 1e-8)
    loss_dict['total_density']['input'].append(input_mass_loss)
    input_str+= 'Total Density Loss: '+'{:0.5f}'.format(input_mass_loss)+'\n'

    pred_mass_loss = abs(act_total_density - pred_total_density) / (act_total_density+ 1e-8)
    loss_dict['total_density']['pred'].append(pred_mass_loss)
    pred_str+= 'Total Density Loss: '+'{:0.5f}'.format(pred_mass_loss)+'\n'

    act_mass_loss = abs(act_total_density - act_total_density) / (act_total_density+ 1e-8)
    act_str+= 'Total Density Loss: '+'{:0.5f}'.format(act_mass_loss)+'\n'

    input_FFTL = F.mse_loss(torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(input_arr)))), torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(act_output_arr)))))
    loss_dict['fft']['input'].append(input_FFTL)
    input_str+= 'FFTL: '+'{:0.5f}'.format(input_FFTL)+'\n'
    pred_FFTL = F.mse_loss(torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(pred_output_arr)))), torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(act_output_arr)))))
    loss_dict['fft']['pred'].append(pred_FFTL)
    pred_str+= 'FFTL: '+'{:0.5f}'.format(pred_FFTL)+'\n'
    actual_FFTL = F.mse_loss(torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(act_output_arr)))), torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(act_output_arr)))))
    act_str+= 'FFTL: '+'{:0.5f}'.format(actual_FFTL)+'\n'

    hd_input = torch.heaviside(torch.from_numpy(input_arr)-torch.quantile(torch.from_numpy(input_arr.flatten()), 0.99), torch.tensor([1.0]))
    hd_pred = torch.heaviside(torch.from_numpy(pred_output_arr)-torch.quantile(torch.from_numpy(pred_output_arr.flatten()), 0.99), torch.tensor([1.0]))
    hd_output = torch.heaviside(torch.from_numpy(act_output_arr)-torch.quantile(torch.from_numpy(act_output_arr.flatten()), 0.99), torch.tensor([1.0]))
    hd_input_loss = F.l1_loss(hd_input, hd_output).sum().item()
    loss_dict['hd']['input'].append(hd_input_loss)
    input_str+= 'HD Loss: '+'{:0.5f}'.format(hd_input_loss)+'\n'
    hd_pred_loss = F.l1_loss(hd_pred, hd_output).sum().item()
    loss_dict['hd']['pred'].append(hd_pred_loss)
    pred_str+= 'HD Loss: '+'{:0.5f}'.format(hd_pred_loss)+'\n'
    actual_HD =F.l1_loss(hd_output, hd_output).sum().item()
    act_str+= 'HD Loss: '+'{:0.5f}'.format(actual_HD)+'\n'

    fig = plt.figure(figsize=(12, 12))
    
    for i in range(3):
        ax  = fig.add_subplot(3, 3, 3*i+1)
        c = ax.pcolormesh(np.mean(input_arr, axis=i).T, vmin=-1, vmax=1, cmap='gnuplot2')
        fig.colorbar(c, ax=ax)
        if i==0:
            ax.set_title('Input: '+str(label[0])+'_'+str(label[1])+'_'+str(label[2]))
        ax.grid('both')
        if i==2:
            ax.text(0.15, -0.02, input_str, fontsize=10, transform=plt.gcf().transFigure, bbox=dict(boxstyle="round", edgecolor='black', facecolor='white'))

    for i in range(3):
        ax2  = fig.add_subplot(3, 3, 3*i+2)
        c = ax2.pcolormesh(np.mean(pred_output_arr, axis=i).T, vmin=-1, vmax=1, cmap='gnuplot2')
        fig.colorbar(c, ax=ax2)
        if i==0:
            ax2.set_title('Predicted Output '+str(label[0])+'_'+str(label[1])+'_'+str(label[2]))
        ax2.grid('both')
        if i==2:
            ax2.text(0.4, -0.02, pred_str, fontsize=10, transform=plt.gcf().transFigure, bbox=dict(boxstyle="round", edgecolor='black', facecolor='white'))

    for i in range(3):
        ax3  = fig.add_subplot(3, 3, 3*i+3)
        c = ax3.pcolormesh(np.mean(act_output_arr, axis=i).T, vmin=-1, vmax=1, cmap='gnuplot2')
        fig.colorbar(c, ax=ax3)
        if i==0:
            ax3.set_title('Actual Output '+str(label[0])+'_'+str(label[1])+'_'+str(label[2]))
        ax3.grid('both')
        if i==2:
            ax3.text(0.7, -0.02, act_str, fontsize=10, transform=plt.gcf().transFigure, bbox=dict(boxstyle="round", edgecolor='black', facecolor='white'))


    if save_filename:
        plt.savefig(save_filename, bbox_inches='tight')
    
    if dont_show==False:
        plt.show()
    plt.close()



def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def set_axis_style(ax, labels):
    ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Sample name')

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


from decimal import Decimal
def compare_output(input_arr, output_arr, labels, args, dir_name ='3D_u502_pred'):
    # Load the trained model
    model = UNet3D(in_ch=1, out_ch=1).to(DEVICE)
    model.load_state_dict(torch.load(args.run_name+'_'+args.loss_type+'_DeeperAttention_'+MODELFILE))
    print(f"Model loaded with {count_parameters(model):,} trainable parameters.", flush=True)
    # Set the model to evaluation mode
    model.eval()

    L1_loss_dict = {'input': [], 'pred': []}
    hist_loss_dict = {'input': [], 'pred': []}
    total_density_loss_dict = {'input': [], 'pred': []}
    fft_loss_dict = {'input': [], 'pred': []}
    hd_loss_dict = {'input': [], 'pred': []}
    Loss_dict = {'l1': L1_loss_dict, 'hist': hist_loss_dict, 'total_density': total_density_loss_dict, 'fft': fft_loss_dict, 'hd': hd_loss_dict}

    with torch.no_grad():
        for i in range(len(input_arr)):
            inputs = input_arr[i].unsqueeze(0).unsqueeze(0)
            targets = output_arr[i].unsqueeze(0).unsqueeze(0)
            label = labels[i]
            inputs, targets = inputs.to(DEVICE), targets.to(DEVICE)
            outputs = model(inputs)
            # Save the output images
            output_image = outputs.cpu().numpy()+ inputs.cpu().numpy()
            plot_comparison(inputs[0][0].numpy(), output_image[0][0], targets[0][0].numpy(), label, Loss_dict, save_filename='./'+dir_name+'/img'+"_".join(str(x) for x in label)+'.png', dont_show = True)
    print("Testing completed and output images plotted.", flush=True)

    pred_data = []
    input_data = []
    Loss_names = Loss_dict.keys()
    for key in Loss_names:
        pred_data.append(np.sort(Loss_dict[key]['pred']))
        input_data.append(np.sort(Loss_dict[key]['input']))
    
    pred_data = np.array(pred_data).T
    input_data = np.array(input_data).T

    fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(24, 16))
    parts = ax1.violinplot(
            pred_data, showmeans=False, showmedians=False,
            showextrema=False)

    for pc in parts['bodies']:
        pc.set_facecolor('plum')
        pc.set_edgecolor('indigo')
        pc.set_alpha(1)

    quartile1, medians, quartile3 = np.percentile(pred_data, [25, 50, 75], axis=0)
    medians_orig= np.percentile(input_data,  50, axis=0)
    whiskers = np.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(pred_data, quartile1, quartile3)])
    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]

    inds = np.arange(1, len(medians) + 1)
    #prediction norm
    ax1.scatter(inds, medians, marker='o', color='white', s=30, zorder=3)
    for i, txt in enumerate(medians):
        ax1.annotate('%.3E'%Decimal(txt), (inds[i]+0.1, medians[i]))
    ax1.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
    ax1.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

    #original as predicted
    ax1.scatter(inds, medians_orig, marker='x', color='red', s=50, zorder=3)
    for i, txt in enumerate(medians_orig):
        ax1.annotate('%.3E'%Decimal(txt), (inds[i]-0.1, medians_orig[i]))
    print(medians_orig)

    ax1.axhline(0, color='black', linestyle='--', lw=2)
    ax1.set_ylabel('Loss value', fontsize=18)
    # set style for the axes
    for ax in [ax1]:
        set_axis_style(ax, Loss_names)
    
    ax1.set_xlabel('Losses', fontsize=18)
    ax1.set_yscale('log')
    ax1.set_ylim([1e-10,2])
    ax1.set_title(args.run_name)


    plt.savefig(args.run_name+'_u502Test.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


parser = argparse.ArgumentParser()
args = parser.parse_args(args=[])
args.batch_size = 8
args.image_size = IMAGESIZE
args.device = DEVICE
args.run_name = "Unet_3D_HD_Delta"
args.epochs = 30
args.lr = 3e-4
args.loss_type = 'dynamic' # 'static' or 'dynamic'


compare_output(input_image_arr, output_image_arr, label_arr, args, dir_name ='3D_u502_pred')

parser = argparse.ArgumentParser()
args = parser.parse_args(args=[])
args.batch_size = 8
args.image_size = IMAGESIZE
args.device = DEVICE
args.run_name = "Unet_3D_HD_Delta_scaled"
args.epochs = 30
args.lr = 3e-4
args.loss_type = 'dynamic' # 'static' or 'dynamic'


compare_output(input_image_arr, output_image_arr, label_arr, args, dir_name ='3D_u502_pred_scaled')
