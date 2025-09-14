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
    for filename in os.listdir(input_dir):
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
    all_files = os.listdir(input_dir)
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
                         num_workers=40,    # adjust based on your CPU
                         pin_memory=False)  # if using GPU

    # Initialize model, loss function, and optimizer
    model = UNet3D(in_ch=1, out_ch=1, args=args).to(DEVICE)
    # For trainable weights
    n_losses = 5  # L1, hist, Mass, Spectral, High Density
    loss_w = torch.nn.Parameter(torch.ones(n_losses, device=DEVICE), requires_grad=True)
    optimizer = optim.AdamW(list(model.parameters()) + [loss_w], lr=args.lr, weight_decay=args.Adamw_weight_decay)
    writer = SummaryWriter(log_dir=f"{args.model_dir}{args.phase}/logs/{args.run_name}")

    start_epoch = 0
    loss_history = []
    #check if model already exists
    if os.path.exists(f"{args.model_dir}{args.phase}/{args.run_name}_{MODELFILE}"):
        checkpoint = torch.load(f"{args.model_dir}{args.phase}/{args.run_name}_{MODELFILE}", map_location=DEVICE, weights_only=False)
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
    for epoch in range(start_epoch, args.epochs):
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
                f"{args.model_dir}{args.phase}/{args.run_name}_{MODELFILE}")
    if args.loss_type != 'static':
        torch.save(loss_w.detach().cpu(), f"{args.model_dir}{args.phase}/plots/loss_data_{args.run_name}_loss.pt")
    import matplotlib.pyplot as plt
    print(loss_history)
    plt.plot(loss_history, marker='o')
    plt.xlabel("Epoch")
    plt.ylabel("Avg Weighted Loss")
    plt.title(f"{args.run_name} Training Loss")
    plt.savefig(f"{args.model_dir}{args.phase}/plots/training_loss_{args.run_name}.png")
    plt.close()

    return min(loss_history)  # return the minimum loss for this run

def model_search_space():
    """
    Define the search space for hyperparameter tuning.
    Returns:
        A dictionary containing the search space parameters.
    """
    search_space_training = {
        'learning_rate': [3e-4, 1e-3, 1e-2],
        'batch_size': [8, 16],
        'Adamw_weight_decay': [1e-2, 1e-4],
        'DWA_temperature': [1, 2]
    }
    search_space_architecture = {
        'base_f': [16, 32, 64],
        'depth': [3, 4, 5],
        'convolution_param': [(3, 1, 1, 1), (3, 1, 2, 2), (3, 1, 3, 3), (4, 1, 3, 2), (5, 1, 2, 1), (5, 1, 4, 2)]        #(kernel size, stride, padding, dilation)
    }

    search_space_tuning = {
        'attention_int_division': [2, 4],
        'LRLU_slope': [0.0, 0.01, 0.1], # LeakyReLU slope, 0 for ReLU
        'dropout_rate': [0.1, 0.2]
    }

    return search_space_training, search_space_architecture, search_space_tuning

def all_models_phase_one():
    """
    Generate all combinations of hyperparameters for phase one, searching across training subspace.
    Returns:
        A table containing a unique combination of hyperparameters as rows.
    """
    search_space_training, search_space_architecture, search_space_tuning = model_search_space()
    from itertools import product

    keys, values = zip(*search_space_training.items())
    combinations = [dict(zip(keys, v)) for v in product(*values)]
    
    import pandas as pd
    df = pd.DataFrame(combinations)
    #Now fix the architecture parameters
    df['model_name'] = df.apply(lambda row: f"model_{row.name + 1}", axis=1)
    df['base_f'] = df.apply(lambda row: 32, axis=1)
    df['depth'] = df.apply(lambda row: 4, axis=1)
    df['attention_int_division'] = df.apply(lambda row: 2, axis=1)
    df['dropout_rate'] = df.apply(lambda row: 0.1, axis=1)
    df['LRLU_slope'] = df.apply(lambda row: 0.0, axis=1)
    df['convolution_param'] = df.apply(lambda row: (3, 1, 1, 1), axis=1)

    return df

def all_models_phase_two():
    """
    Generate all combinations of hyperparameters for phase two, searching across architecture subspace.
    Returns:
        A table containing a unique combination of hyperparameters as rows.
    """
    import pandas as pd
    search_space_training, search_space_architecture, search_space_tuning = model_search_space()
    from itertools import product

    keys, values = zip(*search_space_architecture.items())
    combinations = [dict(zip(keys, v)) for v in product(*values)]
    model_df = pd.read_csv("./phase1_model_performance.csv")
    big_df = pd.DataFrame()
    for i in range(2):
        best_model = model_df.loc[i]
        df = pd.DataFrame(combinations)
        #Now fix the architecture parameters, top 2 phase 1 models
        df['learning_rate'] = df.apply(lambda row: best_model['learning_rate'], axis=1)
        df['batch_size'] = df.apply(lambda row: best_model['batch_size'], axis=1)
        df['Adamw_weight_decay'] = df.apply(lambda row: best_model['Adamw_weight_decay'], axis=1)
        df['DWA_temperature'] = df.apply(lambda row: best_model['DWA_temperature'], axis=1)
        df['attention_int_division'] = df.apply(lambda row: 2, axis=1)
        df['dropout_rate'] = df.apply(lambda row: 0.1, axis=1)
        df['LRLU_slope'] = df.apply(lambda row: 0.0, axis=1)
        #concatenate to big_df
        big_df = pd.concat([big_df, df], ignore_index=True)

    big_df['model_name'] = big_df.apply(lambda row: f"model_{row.name + 1}", axis=1)
    big_df['isTraining'] = big_df.apply(lambda row: 0, axis=1)  #indicating this is model is being currently trained by a running script
    # Add columns for model performance metrics
    big_df['min_training_loss'] = big_df.apply(lambda row: 0, axis=1)
    big_df['avg_test_MSE'] = big_df.apply(lambda row: 0, axis=1)
    big_df['avg_test_density_diff'] = big_df.apply(lambda row: 0, axis=1)
    big_df['training_time_mins'] = big_df.apply(lambda row: 0, axis=1)
    big_df['current_epoch'] = big_df.apply(lambda row: 0, axis=1)
    big_df['keep_model'] = big_df.apply(lambda row: 1, axis=1)

    return big_df

def all_models_phase_three():
    """
    Generate all combinations of hyperparameters for phase three, searching across fine tuning parameters subspace.
    Returns:
        A table containing a unique combination of hyperparameters as rows.
    """
    import pandas as pd
    search_space_training, search_space_architecture, search_space_tuning = model_search_space()
    from itertools import product

    keys, values = zip(*search_space_tuning.items())
    combinations = [dict(zip(keys, v)) for v in product(*values)]
    model_df = pd.read_csv("./phase2_model_performance.csv")
    big_df = pd.DataFrame()
    for i in range(5):
        best_model = model_df.loc[i]
        df = pd.DataFrame(combinations)
        #Now fix the architecture parameters, top 5 phase 2 models
        df['learning_rate'] = df.apply(lambda row: best_model['learning_rate'], axis=1)
        df['batch_size'] = df.apply(lambda row: best_model['batch_size'], axis=1)
        df['Adamw_weight_decay'] = df.apply(lambda row: best_model['Adamw_weight_decay'], axis=1)
        df['DWA_temperature'] = df.apply(lambda row: best_model['DWA_temperature'], axis=1)
        df['base_f'] = df.apply(lambda row: best_model['base_f'], axis=1)
        df['depth'] = df.apply(lambda row: best_model['depth'], axis=1)
        df['convolution_param'] = df.apply(lambda row: best_model['convolution_param'], axis=1)
        #concatenate to big_df
        big_df = pd.concat([big_df, df], ignore_index=True)

    big_df['model_name'] = big_df.apply(lambda row: f"model_{row.name + 1}", axis=1)
    big_df['isTraining'] = big_df.apply(lambda row: 0, axis=1)  #indicating this is model is being currently trained by a running script
    # Add columns for model performance metrics
    big_df['min_training_loss'] = big_df.apply(lambda row: 0, axis=1)
    big_df['avg_test_MSE'] = big_df.apply(lambda row: 0, axis=1)
    big_df['avg_test_density_diff'] = big_df.apply(lambda row: 0, axis=1)
    big_df['training_time_mins'] = big_df.apply(lambda row: 0, axis=1)
    big_df['current_epoch'] = big_df.apply(lambda row: 0, axis=1)
    big_df['keep_model'] = big_df.apply(lambda row: 1, axis=1)

    return big_df