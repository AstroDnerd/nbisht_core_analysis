import os
import copy
import numpy as np
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
from tqdm import tqdm
from torch import optim
from modules import UNet_conditional, EMA
import torchvision
from torchvision.utils import save_image
import logging
from torch.utils.tensorboard import SummaryWriter

logging.basicConfig(format="%(asctime)s - %(levelname)s: %(message)s", level=logging.INFO, datefmt="%I:%M:%S")

class Diffusion:
    def __init__(self, noise_steps=1000, beta_start=1e-4, beta_end=0.02, img_size=64, n_channels=2, device="cuda"):
        self.noise_steps = noise_steps
        self.beta_start = beta_start
        self.beta_end = beta_end
        self.n_channels = n_channels

        self.beta = self.prepare_noise_schedule().to(device)
        self.alpha = 1. - self.beta
        self.alpha_hat = torch.cumprod(self.alpha, dim=0)

        self.img_size = img_size
        self.device = device

    def prepare_noise_schedule(self): #linear scheduler
        return torch.linspace(self.beta_start, self.beta_end, self.noise_steps)

    def noise_images(self, x, t):
        sqrt_alpha_hat = torch.sqrt(self.alpha_hat[t])[:, None, None, None]
        sqrt_one_minus_alpha_hat = torch.sqrt(1 - self.alpha_hat[t])[:, None, None, None]
        Ɛ = torch.randn_like(x)
        return sqrt_alpha_hat * x + sqrt_one_minus_alpha_hat * Ɛ, Ɛ

    def sample_timesteps(self, n):
        return torch.randint(low=1, high=self.noise_steps, size=(n,))

    def sample(self, model, n, initial_condition_image, cfg_scale=3):
        logging.info(f"Sampling {n} new images....")
        model.eval()
        with torch.no_grad():
            xinput = torch.randn((n, int(self.n_channels/2), self.img_size, self.img_size)).to(self.device)
            if cfg_scale > 0:
                x_IC_cfg = torch.normal(mean=0.0, std = 1, size=(n, int(self.n_channels/2), self.img_size, self.img_size)).to(self.device)
            initial_condition_image = torch.reshape(initial_condition_image, xinput.shape)
            x = torch.cat((initial_condition_image,xinput),1)
            for i in tqdm(reversed(range(1, self.noise_steps)), position=0):
                t = (torch.ones(n) * i).long().to(self.device)
                predicted_noise = model(x, t)
                if cfg_scale > 0:
                    xcfg = torch.cat((x_IC_cfg,xinput),1)
                    uncond_predicted_noise = model(xcfg, t)
                    predicted_noise = torch.lerp(uncond_predicted_noise, predicted_noise, cfg_scale)
                alpha = self.alpha[t][:, None, None, None]
                alpha_hat = self.alpha_hat[t][:, None, None, None]
                beta = self.beta[t][:, None, None, None]
                if i > 1:
                    noise = torch.randn_like(x)
                else:
                    noise = torch.zeros_like(x)
                x = 1 / torch.sqrt(alpha) * (x - ((1 - alpha) / (torch.sqrt(1 - alpha_hat))) * predicted_noise) + torch.sqrt(beta) * noise
        model.train()
        x[:,0,:,:] = initial_condition_image
        x = (x.clamp(-1, 1) + 1) / 2
        #x = (x * 255).type(torch.uint8)
        return x

    def get_prediction(self, model, n, initial_condition_dataset, IMAGESIZE = None):
        if IMAGESIZE==None:
            IMAGESIZE = self.img_size
        import numpy as np
        H_initial, edges = np.histogramdd(initial_dataset_core, bins = (IMAGESIZE, IMAGESIZE))
        initial_scalor = max(map(max, H_initial))
        H_initial = (H_initial/initial_scalor)
        sampled_image = self.sample(model, n, H_initial, cfg_scale=0)
        sampled_image[:,0,:,:] = sampled_image[:,0,:,:]*initial_scalor
        sampled_image[:,1,:,:] = sampled_image[:,1,:,:]*initial_scalor*IMAGESIZE
        return sampled_image

        


def setup_logging(run_name):
    os.makedirs("models", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    os.makedirs(os.path.join("models", run_name), exist_ok=True)
    os.makedirs(os.path.join("results", run_name), exist_ok=True)

def save_images(images, edges, path, **kwargs):
    print(path)
    gridlength = images.shape[0]
    fig,ax_arr = plt.subplots(2,gridlength, figsize=(4*gridlength+4,10),  **kwargs)
    for ax_id in range(gridlength):
        if gridlength==1:
            ax0 = ax_arr[0]
            ax1 = ax_arr[1]
        else:
            ax0 = ax_arr[0,ax_id]
            ax1 = ax_arr[1,ax_id]
        z1_plot = ax0.pcolormesh(edges[0], edges[1], images[ax_id][0].T, cmap = 'Grays')
        plt.colorbar(z1_plot,ax=ax0)
        z2_plot = ax1.pcolormesh(edges[0], edges[1], images[ax_id][1].T, cmap = 'Grays')
        if ax_id == 0:
            ax0.set_title("Initial Conditon")
            ax1.set_title("Sink Prediction")
        plt.colorbar(z2_plot,ax=ax1)
    plt.show()
    fig.savefig(path)
    

def train(args, dataloader, edges = None):
    setup_logging(args.run_name)
    device = args.device
    channel = args.channel
    noise_steps = args.noise_steps
    model = UNet_conditional(c_in=channel, c_out=channel, time_dim=256, device=device).to(device)
    optimizer = optim.AdamW(model.parameters(), lr=args.lr)
    mse = nn.MSELoss()
    diffusion = Diffusion(noise_steps = noise_steps, n_channels = channel, img_size=args.image_size, device=device)
    logger = SummaryWriter(os.path.join("runs", args.run_name))
    l = len(dataloader)
    ema = EMA(0.995)
    ema_model = copy.deepcopy(model).eval().requires_grad_(False)

    for epoch in range(args.epochs):
        logging.info(f"Starting epoch {epoch}:")
        pbar = tqdm(dataloader)
        for i, (images, _) in enumerate(pbar):
            images = images.to(device)
            t = diffusion.sample_timesteps(images.shape[0]).to(device)
            x_t, noise = diffusion.noise_images(images, t)
            predicted_noise = model(x_t, t)
            loss = mse(noise, predicted_noise)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            ema.step_ema(ema_model, model)

            pbar.set_postfix(MSE=loss.item())
            logger.add_scalar("MSE", loss.item(), global_step=epoch * l + i)

        if epoch % 10 == 0:
            sampled_images = diffusion.sample(model, initial_condition_image = images[0][0],n=images.shape[0], cfg_scale=0)
            ema_sampled_images = diffusion.sample(ema_model, initial_condition_image = images[0][0],n=images.shape[0], cfg_scale=0)
            save_images(sampled_images, edges, os.path.join("results", args.run_name, f"{epoch}.jpg"))
            save_images(ema_sampled_images, edges, os.path.join("results", args.run_name, f"{epoch}_ema.jpg"))
            torch.save(model.state_dict(), os.path.join("models", args.run_name, f"ckpt.pt"))
            torch.save(ema_model.state_dict(), os.path.join("models", args.run_name, f"ema_ckpt.pt"))
            torch.save(optimizer.state_dict(), os.path.join("models", args.run_name, f"optim.pt"))



