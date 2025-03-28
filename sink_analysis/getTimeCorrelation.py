#import modules
import yt
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import os
import warnings
import random
import sys
warnings.filterwarnings("ignore")
import json
import scipy.fft as fft
from starter2 import *
import track_loader as TL
from animate_plots import animator
import pdb
mp.rcParams.update(mp.rcParamsDefault)
mp.rcParams['agg.path.chunksize'] = 100000
import matplotlib.gridspec as gridspec
import scienceplots

plt.style.use(['science','grid'])

#define paths sinkoff
path_to_data_128 = "/data/cb1/nbisht/anvil_scratch/projects/128/B2/"
path_to_output_plots = "./plots_to_sort/TimeCorrelation/"

#Average Time difference between snapshots
DELTAT = 0.0003744333333

def make_dir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def ac_fft(x):
    N = len(x)
    x_fft = fft.fft(x,n=2*N, axis=0)
    x_fft_conj = np.conj(x_fft)
    autocorr = fft.ifft(x_fft*x_fft_conj, axis=0)
    autocorr = autocorr[:N]
    return np.real(autocorr/autocorr[0])

def ac_self(x):
    N = len(x)
    ac_self = []
    denominator = np.sum(np.multiply(x,x),axis=0)/len(x)
    for i in range(N):
        pdt_arr = np.multiply(x[i:],x[:N-i])
        ac_self.append((np.sum(pdt_arr,axis=0)/len(pdt_arr))/denominator)
    ac_self = np.array(ac_self)
    return ac_self

def get_integral_scale(x):
    L_scale = 0
    for val in x:
        if val>0:
            L_scale+=val
        else:
            return L_scale
    return L_scale

def plot_autocorrelation_0d(x):
    ac_x = ac_fft(x)
    ac_s = ac_self(x)

    fig, (ax1, ax2) = plt.subplots(2, 1)


    ax1.plot(x)
    ax1.set_title('Raw data')

    ax2.plot(ac_x, label='fft')
    ax2.plot(ac_s, label='self')
    ax2.legend()
    ax2.set_title('FFT L:%.3f; Self L:%.3f'%(get_integral_scale(ac_x),get_integral_scale(ac_s)))

    plt.tight_layout()
    plt.savefig("./plots_to_sort/AutocorrelationTest/AutocorrelationTest_0D.png")

def plot_autocorrelation_nd(x, dim='1D'):
    #calculate Nd autocorrelation
    ac_x = ac_fft(x)
    print("FFT shape:",ac_x.shape)
    ac_s = ac_self(x)
    print("self shape:",ac_s.shape)

    ac_x_summed = ac_x
    ac_s_summed = ac_s
    for i in range(1,len(ac_x.shape)):
        print("Compressing along dimension:",i)
        ac_x_summed = np.sum(ac_x_summed, axis=1)/ac_x.shape[1]
        ac_s_summed = np.sum(ac_s_summed, axis=1)/ac_x.shape[1]
    L_fft = get_integral_scale(ac_x_summed)
    L_self = get_integral_scale(ac_s_summed)

    #Make timeseries
    T = np.arange(ac_x.shape[0])

    print(ac_x_summed.shape)
    print(L_fft)
    print(ac_s_summed.shape)
    print(L_self)
    for i in range(len(x)):
        x_i = x[i]
        fig = plt.figure(figsize=(20,20))

        if dim[0] == '1':
            ax1 = fig.add_subplot(211)
            ax1.plot(x_i)
            ax2 = fig.add_subplot(212)
        elif dim[0] =='2':
            ax1 = fig.add_subplot(211)
            ax1.imshow(x_i.T, vmin=-1, vmax=+1)
            ax2 = fig.add_subplot(212)
        elif dim[0]=='3':
            ax1 = fig.add_subplot(221)
            ax1.imshow(np.sum(x_i, axis=0)/x_i.shape[0], vmin=-1, vmax=+1)
            ax1 = fig.add_subplot(222)
            ax1.imshow(np.sum(x_i, axis=1)/x_i.shape[1], vmin=-1, vmax=+1)
            ax1 = fig.add_subplot(223)
            ax1.imshow(np.sum(x_i, axis=2)/x_i.shape[2], vmin=-1, vmax=+1)
            ax2 = fig.add_subplot(224)

        ax1.set_title('Raw data: %04d'%(i))        
        ax2.scatter(T[:i+1], ac_s_summed[:i+1], label='self', c='blue')
        ax2.scatter(T[:i+1], ac_x_summed[:i+1], label='fft', c='orange')
        ax2.set_xlim([0,ac_x_summed.shape[0]])
        ax2.set_ylim([-1,1])
        ax2.legend()

        ax2.set_title('FFT L:%.3f, %.3fs;\n Self L:%.3f, %.3fs'%(L_fft,DELTAT*L_fft,L_self,DELTAT*L_self))

        plt.savefig("./plots_to_sort/AutocorrelationTest/AT_"+dim+"_plots/DD%04d.png"%(i))
        plt.close()
    animator("./plots_to_sort/AutocorrelationTest/AT_"+dim+"_plots","AutocorrelationTest_"+dim)
    
def paperplot_autocorrelation(x):
    #calculate Nd autocorrelation
    ac_x = ac_fft(x)
    print("FFT shape:",ac_x.shape)

    ac_x_summed = ac_x
    for i in range(1,len(ac_x.shape)):
        print("Compressing along dimension:",i)
        ac_x_summed = np.sum(ac_x_summed, axis=1)/ac_x.shape[1]
    L_fft = get_integral_scale(ac_x_summed)

    #Make timeseries
    T = np.arange(ac_x.shape[0])*DELTAT

    print(ac_x_summed.shape)
    print(L_fft)
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111) 
    ax.plot(T, ac_x_summed, c='orange')
    #ax.set_xlim([0,ac_x_summed.shape[0]])
    ax.set_xlabel(r'Time Lag $\tau$ (s)', fontsize=12)
    ax.set_ylabel(r'Autocorrelation Function ACF($\tau$)', fontsize=12)
    ax.hlines(0,0,T[-1],linestyles='dashed', color='black')
    ax.fill_between(T,ac_x_summed, color="none", hatch="X", edgecolor="powderblue", linewidth=0.0)
    ax.set_ylim([-0.1,1])
    props = dict(boxstyle='round', facecolor='wheat', alpha = 0.6)
    # place a text box in upper left in axes coords
    ax.text(0.15, 0.2, 'L: %.3fs'%(DELTAT*L_fft), transform=ax.transAxes, fontsize=14,
            verticalalignment='bottom', bbox=props, zorder=5)
    print('FFT L:%.3f, %.3fs'%(L_fft,DELTAT*L_fft))

    plt.savefig("./plots_to_sort/Densityvar_autocorrelation.png")
    plt.close()

#Test Cases
if 0:
    #0D, 1 number varying through time from -1 to 1, sinusoidal at first then random
    x = np.linspace(0, 2 * np.pi, 50)
    x = np.sin(x)
    x = np.append(x,np.random.uniform(-1, 1, 75))
    plot_autocorrelation_0d(x)
if 0:
    #1D, initial a 1D sine curve then randomness
    x = []
    for i in range(125):
        if i<50:
            shift_i = i*np.pi/25
            x_i = np.linspace(0+shift_i, 2 * np.pi+shift_i, 128)
            x.append(np.sin(x_i))
        else:
            x.append(np.random.uniform(-1, 1, 128))
    plot_autocorrelation_nd(x)
if 0:
    #2D, initial a 2D sine sheet then randomness
    x = []
    for i in range(125):
        if i<50:
            shift_i = i*np.pi/25
            x_i = np.linspace(0+shift_i, 2 * np.pi+shift_i, 128)
            x_i = np.tile(x_i,(128,1))
            x.append(np.sin(x_i))
        else:
            x_i = np.random.uniform(-1, 1, (128,128))
            x.append(x_i)
    plot_autocorrelation_nd(x, dim='2D')
if 0:
    #3D, initial a 3D sine sheet then randomness
    x = []
    for i in range(125):
        if i<50:
            shift_i = i*np.pi/25
            x_i = np.linspace(0+shift_i, 2 * np.pi+shift_i, 128)
            x_i = np.repeat(np.tile(x_i,(128,1))[:,:,np.newaxis],128,axis=2)
            x.append(np.sin(x_i))
        else:
            x_i = np.random.uniform(-1, 1, (128,128,128))
            x.append(x_i)
    plot_autocorrelation_nd(x, dim='3D')
if 1:
    #3D density Timeseries_Cube
    nonsink_trackname = 'nb101'
    this_track = track_info.tracks[nonsink_trackname]
    df_name = this_track.sim_directory+'/datasets/nb101_TimeseriesCubes_Density.npy'
    infile = open(df_name, 'rb')
    TSCube_density = np.load(infile)
    time_array = np.load(infile)
    infile.close()
    #density variance
    TSCube_density = TSCube_density-1
    #plot_autocorrelation_nd(TSCube_density, dim='3D_densityvariance')
    paperplot_autocorrelation(TSCube_density)

if 0:
    #3D density Timeseries_Cube visualised with plotly
    import plotly.graph_objects as go
    import plotly.express as px
    nonsink_trackname = 'nb101'
    this_track = track_info.tracks[nonsink_trackname]
    df_name = this_track.sim_directory+'/datasets/nb101_TimeseriesCubes_Density.npy'
    infile = open(df_name, 'rb')
    TSCube_density = np.load(infile)
    time_array = np.load(infile)
    infile.close()
    html_name = this_track.sim_directory+'/datasets/nb101_TimeseriesCubes_Density_visualiser_DD0125.html'


    #plotly figure, first create x,y,z meshgrid with evenly spaced numbers
    X, Y, Z = np.mgrid[0:1:128j, 0:1:128j, 0:1:128j]
    #flatten 3D array
    flat_tscube = TSCube_density[120].flatten()
    minval = np.quantile(flat_tscube,0.2)
    base_minval = 10**(round(np.log10(minval)))
    isominval = int(minval/base_minval)*base_minval
    maxval = np.quantile(flat_tscube,0.99999)
    base_maxval = 10**(round(np.log10(maxval)))
    isomaxval = int(maxval/base_maxval)*base_maxval
    print(minval,maxval)
    print("min isocontour:%0.3f max isocontour:%0.3f"%(isominval,isomaxval))
    fig = go.Figure(data=go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=flat_tscube,
        isomin=isominval,
        isomax=isomaxval,
        opacity=0.1, # needs to be small to see through all surfaces
        surface_count=50, # needs to be a large number for good volume rendering
        ))
    fig.write_html(html_name)


