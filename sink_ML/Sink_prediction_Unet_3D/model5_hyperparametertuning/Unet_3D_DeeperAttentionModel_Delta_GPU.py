# standard system modules
import os, sys
import h5py 
import argparse
# standard module for tabular data
import pandas as pd

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
from model4_modules import *

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
def train_unet(input_arr, output_arr, labels, args):
    X = torch.stack(input_arr, dim=0).float().unsqueeze(1)  # [N,1,64, 64, 64]
    Y = torch.stack(output_arr, dim=0).float().unsqueeze(1) # [N,1,64, 64, 64]

    dataset = TensorDataset(X, Y)
    loader  = DataLoader(dataset,
                         batch_size=args.batch_size,
                         shuffle=True,
                         num_workers=12,    # adjust based on your CPU
                         pin_memory=False)  # if using GPU

    # Initialize model, loss function, and optimizer
    model = UNet3D(in_ch=1, out_ch=1).to(DEVICE)
    # For trainable weights
    n_losses = 5  # L1, hist, Mass, Spectral, High Density
    loss_w = torch.nn.Parameter(torch.ones(n_losses, device=DEVICE), requires_grad=True)
    optimizer = optim.Adam(list(model.parameters()) + [loss_w], lr=args.lr)
    writer = SummaryWriter(log_dir=f"./logs/{args.run_name}")
    L1_scale = 5.0  # scale L1 loss to be more influential
    hist_scale = 10.0  # scale histogram loss
    mass_scale = 0.1   # give mass only 10% initial influence
    spectral_scale = 1.0  # scale spectral loss
    hd_scale = 10.0 #high density loss scale

    #NOW TRAINING ON THE DIFFERENCE RATHER THAN DIRECTLY THE OUTPUT

    if args.loss_type == 'static':
        #Precompute initial-loss weights
        with torch.no_grad():
            x0 = X[0:1].to(DEVICE)
            y0 = Y[0:1].to(DEVICE)
            delta0 = y0-x0
            deltap0 = model(x0)
            #L1
            L1_0   = F.l1_loss(deltap0, delta0).item()
            # Histogram loss
            hist_0 = hist_loss(deltap0, delta0, bins=32)
            #Mass
            mass_true = delta0.sum().item()
            mass_pred = deltap0.sum().item()
            MASS_0 = abs(mass_pred - mass_true)
            MASS_0 = torch.log1p(MASS_0)  # log1p to avoid large values
            #FFT Loss
            spectral_0 = compute_spectral_loss(deltap0, delta0).item()
            #High desnity
            hd_true0 = torch.heaviside(delta0-torch.quantile(delta0.flatten(), 0.99), torch.tensor([1.0]))
            hd_pred0 = torch.heaviside(deltap0-torch.quantile(deltap0.flatten(), 0.99), torch.tensor([1.0]))
            hd_loss0 = F.l1_loss(hd_pred0, hd_true0)
            
            init_w = torch.tensor([1/L1_0, 1/hist_0, 1/MASS_0, 1/spectral_0, 1/hd_loss0],device=DEVICE)
            static_w = init_w / init_w.sum()
    if args.loss_type == 'dynamic':
        prev_losses = [1.0]*n_losses    # dummy for epoch 0
        prev2_losses = [1.0]*n_losses   # dummy for epoch -1
        w = [1.0]*n_losses

    # === Training loop ===
    loss_history = []
    for epoch in range(args.epochs):
        model.train()
        total_epoch_loss = 0.0
        sum_l1, sum_hist, sum_mass, sum_spectral, sum_hd = 0.0, 0.0, 0.0, 0.0, 0.0

        for x, y in loader:
            x, y = x.to(DEVICE), y.to(DEVICE)
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
                T = 2.0 #defined as temperature in the paper, smoothes extremes
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
              f"L1={avg_l1:.4f} Hist={avg_hist:.4f} "
              f"Mass={avg_mass:.4f} "
              f"Spectral={avg_spectral:.4f} "
              f"HighDensity={avg_hd:.4f} "
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
    torch.save(model.state_dict(), args.run_name+'_'+args.loss_type+'_DeeperAttention_'+MODELFILE)
    if args.loss_type != 'static':
        torch.save(loss_w.detach().cpu(), "./plots/deltaloss_weights_deeper_attention.pt")
    import matplotlib.pyplot as plt
    print(loss_history)
    plt.plot(loss_history, marker='o')
    plt.xlabel("Epoch")
    plt.ylabel("Avg Weighted Loss")
    plt.title(f"{args.run_name} Training Loss")
    plt.savefig(f"./plots/{args.run_name}_deeper_attention_training_loss.png")
    plt.close()

#launch training
parser = argparse.ArgumentParser()
args = parser.parse_args(args=[])
args.batch_size = 8
args.image_size = IMAGESIZE
args.device = DEVICE
args.run_name = "Unet_3D_HD_Delta_scaled"
args.epochs = 30
args.lr = 3e-4
args.loss_type = 'dynamic' # 'static' or 'dynamic'

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

    input_FFTL = F.mse_loss(torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(input_arr)))), torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(act_output_arr)))))
    input_str+= 'FFTL: '+'{:0.5f}'.format(input_FFTL)+'\n'
    pred_FFTL = F.mse_loss(torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(pred_output_arr)))), torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(act_output_arr)))))
    pred_str+= 'FFTL: '+'{:0.5f}'.format(pred_FFTL)+'\n'
    actual_FFTL = F.mse_loss(torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(act_output_arr)))), torch.log1p(torch.abs(torch.fft.fftn(torch.from_numpy(act_output_arr)))))
    act_str+= 'FFTL: '+'{:0.5f}'.format(actual_FFTL)+'\n'

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

def compare_output(input_arr, output_arr, labels, args):
    # Load the trained model
    model = UNet3D(in_ch=1, out_ch=1).to(DEVICE)
    model.load_state_dict(torch.load(args.run_name+'_'+args.loss_type+'_DeeperAttention_'+MODELFILE))
    print(f"Model loaded with {count_parameters(model):,} trainable parameters.", flush=True)
    # Set the model to evaluation mode
    model.eval()

    with torch.no_grad():
        for i in range(len(input_arr)):
            inputs = input_arr[i].unsqueeze(0).unsqueeze(0)
            targets = output_arr[i].unsqueeze(0).unsqueeze(0)
            label = labels[i]
            inputs, targets = inputs.to(DEVICE), targets.to(DEVICE)
            outputs = model(inputs)
            # Save the output images
            output_image = outputs.cpu().numpy()+ inputs.cpu().numpy()
            plot_comparison(inputs[0][0].numpy(), output_image[0][0], targets[0][0].numpy(), label, save_filename='./3D_test_plots_scaled/img'+"_".join(str(x) for x in label)+'.png', dont_show = True)
    print("Testing completed and output images plotted.", flush=True)

compare_output(input_image_test, output_image_test, label_test, args)