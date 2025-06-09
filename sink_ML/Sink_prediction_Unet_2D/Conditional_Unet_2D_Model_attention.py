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
from modules import *

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

print(f'Available device: {str(DEVICE):4s}', flush=True)


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
input_image_train, input_image_test, output_image_train, output_image_test, label_train, label_test = train_test_split(input_image_arr, output_image_arr, label_arr, 
                                                                                                                       shuffle = False, test_size=TEST_PERCENTAGE, random_state=seed)
    

import torch.nn as nn
from torchmetrics.image import StructuralSimilarityIndexMeasure
from torch.utils.tensorboard import SummaryWriter

def compute_spectral_loss(pred, target):
    #pred, target are shaped [B, 1, H, W]
    pred_fft = torch.fft.fft2(pred.squeeze(1))   # [B, H, W], complex
    targ_fft = torch.fft.fft2(target.squeeze(1))
    #find the absolute magnitude, whihc gives the power spectrum
    pred_mag = torch.abs(pred_fft)
    targ_mag = torch.abs(targ_fft)
    # power spectrum difference (normalized optional)
    loss = F.mse_loss(pred_mag, targ_mag)
    return loss


class UNet2D(nn.Module):
    def __init__(self, in_channels, out_channels, base_filters=32):
        super().__init__()
        # Encoder path
        self.enc1 = nn.Sequential(
            nn.Conv2d(in_channels, base_filters, kernel_size=3, padding=1),
            nn.ReLU(inplace=True),
            nn.Conv2d(base_filters, base_filters, kernel_size=3, padding=1),
            nn.ReLU(inplace=True),
        )
        self.pool = nn.MaxPool2d(kernel_size=2, stride=2)
        self.enc2 = nn.Sequential(
            nn.Conv2d(base_filters, base_filters*2, kernel_size=3, padding=1),
            nn.ReLU(inplace=True),
            nn.Conv2d(base_filters*2, base_filters*2, kernel_size=3, padding=1),
            nn.ReLU(inplace=True),
        )
        # (You can add more encoder levels here…)

        # Decoder path
        self.upconv2 = nn.ConvTranspose2d(base_filters*2, base_filters, kernel_size=2, stride=2)
        self.dec2 = nn.Sequential(
            nn.Conv2d(base_filters*2, base_filters, kernel_size=3, padding=1),
            nn.ReLU(inplace=True),
            nn.Conv2d(base_filters, base_filters, kernel_size=3, padding=1),
            nn.ReLU(inplace=True),
        )
        # (Add more decoder levels to mirror the encoder…)

        # Final 1×1 conv to map to desired output channels
        self.outconv = nn.Conv2d(base_filters, out_channels, kernel_size=1)

    def forward(self, x):
        # Encoder
        e1 = self.enc1(x)                          # [B, F, H, W]
        e2 = self.enc2(self.pool(e1))              # [B, 2F, H/2, W/2]
        # (encode deeper levels…)

        # Decoder
        d2 = self.upconv2(e2)                      # [B, F, H, W]
        d2 = torch.cat([d2, e1], dim=1)            # skip connection
        d2 = self.dec2(d2)                         # [B, F, H, W]
        # (decode through further levels…)

        return self.outconv(d2)                    # [B, out_channels, H, W]



def train_unet(input_arr, output_arr, labels, args):
    X = torch.stack(input_arr, dim=0).float().unsqueeze(1)  # [N,1,128,128]
    Y = torch.stack(output_arr, dim=0).float().unsqueeze(1) # [N,1,128,128]

    dataset = TensorDataset(X, Y)
    loader  = DataLoader(dataset,
                         batch_size=args.batch_size,
                         shuffle=True,
                         num_workers=8,    # adjust based on your CPU
                         pin_memory=False)  # if using GPU

    # Initialize model, loss function, and optimizer
    model = UNet2D(in_channels=1, out_channels=1).to(DEVICE)
    # For trainable weights
    loss_w = torch.nn.Parameter(torch.ones(5, device=DEVICE), requires_grad=True)
    optimizer = optim.Adam(list(model.parameters()) + [loss_w], lr=args.lr)
    writer = SummaryWriter(log_dir=f"./logs/{args.run_name}")
    mass_scale = 0.1   # give mass only 10% initial influence

    if args.loss_type == 'static':
        #Precompute initial-loss weights
        with torch.no_grad():
            x0 = X[0:1].to(DEVICE)
            y0 = Y[0:1].to(DEVICE)
            p0 = model(x0)
            #L1
            L1_0   = F.l1_loss(p0, y0).item()
            # MSE
            mse = F.mse_loss(p0, y0)
            #SSIM (1 - SSIM)
            ssim_func = StructuralSimilarityIndexMeasure()
            SSIM_0 = 1 - ssim_func(p0, y0)
            #Mass (invert log10 before summation)
            mass_true = torch.pow(10, y0).sum().item()
            mass_pred = torch.pow(10, p0).sum().item()
            MASS_0 = abs(mass_pred - mass_true) / (mass_true+ 1e-8)
            #FFT Loss
            spectral_0 = compute_spectral_loss(p0, y0).item()
            
            init_w = torch.tensor([1/L1_0, 1/mse, 1/SSIM_0, mass_scale/MASS_0, 1/spectral_0],device=DEVICE)
            static_w = init_w / init_w.sum()

    # === Training loop ===
    loss_history = []
    for epoch in range(args.epochs):
        model.train()
        total_epoch_loss = 0.0
        sum_l1, sum_mse, sum_ssim, sum_mass, sum_spectral = 0.0, 0.0, 0.0, 0.0, 0.0

        for x, y in loader:
            x, y = x.to(DEVICE), y.to(DEVICE)
            optimizer.zero_grad()
            pred = model(x)
            #Computing individual losses
            #L1 loss
            l1 = F.l1_loss(pred, y)
            # MSE
            mse = F.mse_loss(pred, y)
            #SSIM loss = 1 − SSIM
            ssim_func = StructuralSimilarityIndexMeasure()
            ssim_l = 1 - ssim_func(pred, y)
            #Mass loss (invert log10)
            mass_t = torch.pow(10, y).sum()
            mass_p = torch.pow(10, pred).sum()
            mass_error = torch.abs(mass_p - mass_t) / (mass_t + 1e-8)
            mass_l = torch.log1p(mass_error)  # log1p to avoid large values
            # Spectral loss using power spectrum
            spectral_l = compute_spectral_loss(pred, y)

            #Stacking losses
            losses = torch.stack([l1, mse, ssim_l, mass_l, spectral_l])
            if args.loss_type == 'static':
                total_loss = (static_w * losses).sum()
            else:
                w = torch.softmax(loss_w, dim=0)
                total_loss = (w * losses).sum()

            total_loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)  # optional
            optimizer.step()
            # accumulate
            batch_size = x.size(0)
            total_epoch_loss += total_loss.item() * batch_size
            sum_l1     += l1.item()      * batch_size
            sum_mse    += mse.item()     * batch_size
            sum_ssim   += ssim_l.item()  * batch_size
            sum_mass   += mass_l.item()  * batch_size
            sum_spectral   += spectral_l.item()  * batch_size

        # compute per-sample averages
        N = len(dataset)
        avg_total = total_epoch_loss / N
        avg_l1    = sum_l1     / N
        avg_mse   = sum_mse    / N
        avg_ssim  = sum_ssim   / N
        avg_mass  = sum_mass   / N
        avg_spectral = sum_spectral / N

        # log to console
        print(f"Epoch {epoch+1}: "
              f"Total={avg_total:.4f} | "
              f"L1={avg_l1:.4f} MSE={avg_mse:.4f} "
              f"SSIM={avg_ssim:.4f} Mass={avg_mass:.4f}"
              f"Spectral={avg_spectral:.4f}")

        # log to TensorBoard
        writer.add_scalar("Loss/Total",    avg_total, epoch)
        writer.add_scalar("Loss/L1",       avg_l1,    epoch)
        writer.add_scalar("Loss/MSE",     avg_mse,   epoch)
        writer.add_scalar("Loss/SSIM",     avg_ssim,  epoch)
        writer.add_scalar("Loss/Mass",     avg_mass,  epoch)
        writer.add_scalar("Loss/Spectral", avg_spectral, epoch)

        if args.loss_type!='static':
            # log the learned weights (after softmax)
            curr_w = torch.softmax(loss_w, dim=0).detach().cpu().tolist()
            for i, name in enumerate(["L1","MSE","SSIM","Mass", "Spectral"]):
                writer.add_scalar(f"Weights/{name}", curr_w[i], epoch)

    writer.close()

    # === Save & plot ===
    torch.save(model.state_dict(), MODELFILE)
    if args.loss_type != 'static':
        torch.save(loss_w.detach().cpu(), "./plots/loss_weights.pt")
    import matplotlib.pyplot as plt
    plt.plot(loss_history, marker='o')
    plt.xlabel("Epoch")
    plt.ylabel("Avg Weighted Loss")
    plt.title(f"{args.run_name} Training Loss")
    plt.savefig(f"./plots/{args.run_name}_training_loss.png")
    plt.close()

#launch training
parser = argparse.ArgumentParser()
args = parser.parse_args(args=[])
args.batch_size = 1
args.image_size = IMAGESIZE
args.device = DEVICE
args.run_name = "Unet_2D_spectral"
args.epochs = 30
args.lr = 3e-4
args.loss_type = 'static' # 'static' or 'dynamic'

print("Starting training!", flush=True)
train_unet(input_image_train, output_image_train, label_train, args)

def plot_comparison(input_arr, pred_output_arr, act_output_arr, label, save_filename='xy_data.png', dont_show = False):
    input_str = 'INPUT IMAGE:\n'
    pred_str = 'PREDICTED IMAGE:\n'
    act_str = 'ACTUAL IMAGE:\n'
    # Calculate metrics
    input_MSE = ski.metrics.mean_squared_error(act_output_arr, input_arr)
    input_str+= 'MSE: '+'{:0.5f}'.format(input_MSE)+'\n'
    pred_MSE = ski.metrics.mean_squared_error(act_output_arr, pred_output_arr)
    pred_str+= 'MSE: '+'{:0.5f}'.format(pred_MSE)+'\n'
    actual_MSE = ski.metrics.mean_squared_error(act_output_arr, act_output_arr)
    act_str+= 'MSE: '+'{:0.5f}'.format(actual_MSE)+'\n'

    input_NMI = ski.metrics.normalized_mutual_information(act_output_arr, input_arr, bins=100)
    input_str+= 'NMI: '+'{:0.5f}'.format(input_NMI)+'\n'
    pred_NMI = ski.metrics.normalized_mutual_information(act_output_arr, pred_output_arr, bins=100)
    pred_str+= 'NMI: '+'{:0.5f}'.format(pred_NMI)+'\n'
    actual_NMI = ski.metrics.normalized_mutual_information(act_output_arr, act_output_arr, bins=100)
    act_str+= 'NMI: '+'{:0.5f}'.format(actual_NMI)+'\n'

    inversed_input_arr = img_inverse_transform(input_arr)
    inversed_pred_output_arr = img_inverse_transform(pred_output_arr)
    inversed_act_output_arr = img_inverse_transform(act_output_arr)
    input_PSNR = ski.metrics.peak_signal_noise_ratio(inversed_act_output_arr, inversed_input_arr, data_range=np.max(inversed_input_arr) - np.min(inversed_input_arr))
    input_str+= 'PSNR: '+'{:0.5f}'.format(input_PSNR)+'\n'
    pred_PSNR = ski.metrics.peak_signal_noise_ratio(inversed_act_output_arr, inversed_pred_output_arr, data_range=np.max(inversed_pred_output_arr) - np.min(inversed_pred_output_arr))
    pred_str+= 'PSNR: '+'{:0.5f}'.format(pred_PSNR)+'\n'
    actual_PSNR = ski.metrics.peak_signal_noise_ratio(inversed_act_output_arr, inversed_act_output_arr, data_range=np.max(inversed_act_output_arr) - np.min(inversed_act_output_arr))
    act_str+= 'PSNR: '+'{:0.5f}'.format(actual_PSNR)+'\n'

    input_SSI = ski.metrics.structural_similarity(act_output_arr, input_arr, gradient=False, data_range=np.max(input_arr) - np.min(input_arr), channel_axis=None)
    input_str+= 'SSI: '+'{:0.5f}'.format(input_SSI)+'\n'
    pred_SSI = ski.metrics.structural_similarity(act_output_arr, pred_output_arr, gradient=False, data_range=np.max(pred_output_arr) - np.min(pred_output_arr), channel_axis=None)
    pred_str+= 'SSI: '+'{:0.5f}'.format(pred_SSI)+'\n'
    actual_SSI = ski.metrics.structural_similarity(act_output_arr, act_output_arr, gradient=False, data_range=np.max(act_output_arr) - np.min(act_output_arr), channel_axis=None)
    act_str+= 'SSI: '+'{:0.5f}'.format(actual_SSI)+'\n'

    input_total_density = np.sum(inversed_input_arr)
    input_str+= 'Total Density: '+'{:0.5f}'.format(input_total_density)+'\n'
    pred_total_density = np.sum(inversed_pred_output_arr)
    pred_str+= 'Total Density: '+'{:0.5f}'.format(pred_total_density)+'\n'
    act_total_density = np.sum(inversed_act_output_arr)
    act_str+= 'Total Density: '+'{:0.5f}'.format(act_total_density)+'\n'

    input_FFTL = F.mse_loss(torch.abs(torch.fft.fft2(input_arr)), torch.abs(torch.fft.fft2(act_output_arr)))
    input_str+= 'FFTL: '+'{:0.5f}'.format(input_FFTL)+'\n'
    pred_FFTL = F.mse_loss(torch.abs(torch.fft.fft2(pred_output_arr)), torch.abs(torch.fft.fft2(act_output_arr)))
    pred_str+= 'FFTL: '+'{:0.5f}'.format(pred_FFTL)+'\n'
    actual_FFTL = F.mse_loss(torch.abs(torch.fft.fft2(act_output_arr)), torch.abs(torch.fft.fft2(act_output_arr)))
    act_str+= 'FFTL: '+'{:0.5f}'.format(actual_FFTL)+'\n'

    fig = plt.figure(figsize=(12, 12))
    
    for i in range(3):
        ax  = fig.add_subplot(3, 3, 3*i+1)
        c = ax.pcolormesh(np.mean(input_arr, axis=i).T)
        fig.colorbar(c, ax=ax)
        if i==0:
            ax.set_title('Input: '+str(label[0])+'_'+str(label[1])+'_'+str(label[2]))
        ax.grid('both')
        if i==2:
            ax.text(0.15, -0.02, input_str, fontsize=10, transform=plt.gcf().transFigure, bbox=dict(boxstyle="round", edgecolor='black', facecolor='white'))

    for i in range(3):
        ax2  = fig.add_subplot(3, 3, 3*i+2)
        c = ax2.pcolormesh(np.mean(pred_output_arr, axis=i).T)
        fig.colorbar(c, ax=ax2)
        if i==0:
            ax2.set_title('Predicted Output '+str(label[0])+'_'+str(label[1])+'_'+str(label[2]))
        ax2.grid('both')
        if i==2:
            ax2.text(0.4, -0.02, pred_str, fontsize=10, transform=plt.gcf().transFigure, bbox=dict(boxstyle="round", edgecolor='black', facecolor='white'))

    for i in range(3):
        ax3  = fig.add_subplot(3, 3, 3*i+3)
        c = ax3.pcolormesh(np.mean(act_output_arr, axis=i).T)
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

def compare_output(input_arr, output_arr, labels, args):
    # Load the trained model
    model = UNet3D(in_channels=1, out_channels=1).to(DEVICE)
    model.load_state_dict(torch.load(args.run_name+'_'+args.loss_type+'_'+MODELFILE))
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
            plot_comparison(inputs[0][0].numpy(), output_image[0][0], targets[0][0].numpy(), label, save_filename='./3D_test_plots/img'+"_".join(str(x) for x in label)+'.png', dont_show = True)
    print("Testing completed and output images plotted.", flush=True)

compare_output(input_image_test, output_image_test, label_test, args)