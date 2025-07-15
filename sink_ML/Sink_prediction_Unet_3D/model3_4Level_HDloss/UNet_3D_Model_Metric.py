# standard system modules
import os, sys
import h5py 
import argparse
# standard module for tabular data
import json

# standard module for array manipulation
import numpy as np
from itertools import permutations

# standard statistical module
import scipy.stats as st
from scipy import linalg
from scipy.stats import ks_2samp
import skimage as ski
from skimage.metrics import structural_similarity as ski_ssim, normalized_mutual_information as ski_nmi

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

DATAFILE  = '/home/x-nbisht1/scratch/projects/128/B2/datasets/nb101_ML_dataset_AllData_AutoEnc.h5'
IMAGEINPUT = '/home/x-nbisht1/scratch/projects/128/B2/datasets/Unet/Images_Input_3D/'
IMAGEOUTPUT = '/home/x-nbisht1/scratch/projects/128/B2/datasets/Unet/Images_Output_3D/'
CORESET  = '/home/x-nbisht1/scratch/projects/128/B2/datasets/nb101_all_frames.h5'
DENSITYCUBES = '/home/x-nbisht1/scratch/projects/128/B2/datasets/nb101_TimeseriesCubes_Density.npy'
MODELFILE = 'nnmodel.dict'

IMAGESIZE = 64
FRAMES = np.arange(0,90, 1)
FRAME_DIFF = 30
TEST_PERCENTAGE = 0.01

DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
#DEVICE = torch.device('cpu')

print(f'Available device: {str(DEVICE):4s}')

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

import torch.nn as nn
from torch.utils.tensorboard import SummaryWriter

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


input_image_dic = get_images(IMAGEINPUT)
output_image_dic = get_images(IMAGEOUTPUT)
input_image_arr = []
output_image_arr = []
label_arr = []
for key in input_image_dic.keys():
    input_image_arr.append(torch.tensor(input_image_dic[key]['img']))
    output_image_arr.append(torch.tensor(output_image_dic[key]['img']))
    label_arr.append(input_image_dic[key]['label'])
input_image_train, input_image_test, output_image_train, output_image_test, label_train, label_test = train_test_split(input_image_arr, output_image_arr, label_arr, 
                                                                                                                       shuffle = False, test_size=TEST_PERCENTAGE, random_state=seed)

def plot_comparison(input_arr, pred_output_arr, act_output_arr, label, save_filename='xy_data.png', dont_show = False):
    # Calculate metrics


    input_flat = np.ravel(input_arr)
    pred_flat = np.ravel(pred_output_arr) 
    actual_flat = np.ravel(act_output_arr)


    fig = plt.figure(figsize=(12, 8))
    ax  = fig.add_subplot(3, 1, 1)
    counts, bins, patches = ax.hist(input_flat, bins=16, alpha=0.5, range = (-4, 4))
    for count, bin_edge, patch in zip(counts, bins, patches):
        # Get the x-coordinate (center of the bar) and y-coordinate (top of the bar)
        x = patch.get_x() + patch.get_width() / 2
        y = patch.get_height()

        # Add the text label slightly above the bar
        if count > 0:  # Only display for non-empty bins
            ax.text(x, y + 0.05 * y, int(count), ha='center', va='bottom', fontsize=8)

    ax.set_title('Input: '+str(label[0])+'_'+str(label[1])+'_'+str(label[2]))
    ax.grid('both')


    ax2  = fig.add_subplot(3, 1, 2)
    counts, bins, patches = ax2.hist(pred_flat, bins=16, alpha=0.5, range = (-4, 4))
    ax2.set_title('Predicted Output '+str(label[0])+'_'+str(label[1])+'_'+str(label[2]))
    for count, bin_edge, patch in zip(counts, bins, patches):
        # Get the x-coordinate (center of the bar) and y-coordinate (top of the bar)
        x = patch.get_x() + patch.get_width() / 2
        y = patch.get_height()

        # Add the text label slightly above the bar
        if count > 0:  # Only display for non-empty bins
            ax2.text(x, y + 0.05 * y, int(count), ha='center', va='bottom', fontsize=8)

    ax2.grid('both')

    ax3  = fig.add_subplot(3, 1, 3)
    counts, bins, patches = ax3.hist(actual_flat, bins=16, alpha=0.5, range = (-4, 4))
    for count, bin_edge, patch in zip(counts, bins, patches):
        # Get the x-coordinate (center of the bar) and y-coordinate (top of the bar)
        x = patch.get_x() + patch.get_width() / 2
        y = patch.get_height()

        # Add the text label slightly above the bar
        if count > 0:  # Only display for non-empty bins
            ax3.text(x, y + 0.05 * y, int(count), ha='center', va='bottom', fontsize=8)

    ax3.set_title('Actual Output '+str(label[0])+'_'+str(label[1])+'_'+str(label[2]))
    ax3.grid('both')

    if save_filename:
        plt.savefig(save_filename, bbox_inches='tight')
    
    if dont_show==False:
        plt.show()
    plt.close()

def compare_output(input_arr, output_arr, labels, args, dir_name ='3D_model_perf'):
    # Load the trained model
    model = UNet3D(in_ch=1, out_ch=1).to(DEVICE)
    model.load_state_dict(torch.load(args.run_name+'_'+args.loss_type+'_DeeperAttention_'+MODELFILE))
    print(f"Model loaded with {count_parameters(model):,} trainable parameters.")
    model.eval()

    with torch.no_grad():
        for i in range(len(input_arr)):
            inputs = input_arr[i].unsqueeze(0).unsqueeze(0)
            targets = output_arr[i].unsqueeze(0).unsqueeze(0)
            label = labels[i]
            inputs, targets = inputs.to(DEVICE), targets.to(DEVICE)
            outputs = model(inputs)
            # Save the output images
            output_image = outputs.cpu().numpy()
            plot_comparison(inputs[0][0].numpy(), output_image[0][0], targets[0][0].numpy(), label, save_filename='./'+dir_name+'/img'+"_".join(str(x) for x in label)+'.png', dont_show = True)
    print("Testing completed and output images plotted.")

parser = argparse.ArgumentParser()
args = parser.parse_args(args=[])
args.batch_size = 8
args.image_size = IMAGESIZE
args.device = DEVICE
args.run_name = "Unet_3D_HD"
args.epochs = 30
args.lr = 3e-4
args.loss_type = 'dynamic' # 'static' or 'dynamic'


compare_output(input_image_test, output_image_test, label_test, args, dir_name ='3D_model_perf')
