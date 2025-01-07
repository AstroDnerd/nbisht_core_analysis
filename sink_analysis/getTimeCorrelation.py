#import modules
import yt
import matplotlib as mpl
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
from animate_plots import animator
import pdb

#define paths sinkoff
path_to_data_128 = "/data/cb1/nbisht/anvil_scratch/projects/128/B2/"
path_to_output_plots = "./plots_to_sort/TimeCorrelation/"

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

def plot_autocorrelation_nd(x):
    #calculate Nd autocorrelation
    ac_x = ac_fft(x)
    print("FFT shape:",ac_x.shape)
    ac_s = ac_self(x)
    print("self shape:",ac_s.shape)

    #Make timeseries
    t = np.arange(ac_x.shape[0])
    y = np.arange(ac_x.shape[1])
    Y,T = np.meshgrid(y,t)
    print(T.shape,Y.shape)
    L_fft = get_integral_scale(np.sum(ac_x, axis=1)/ac_x.shape[1])
    L_self = get_integral_scale(np.sum(ac_s, axis=1)/ac_s.shape[1])
    print(ac_x)
    for i in range(len(x)):
        x_i = x[i]
        fig = plt.figure(figsize=(10,20))

        ax1 = fig.add_subplot(211)
        ax1.plot(x_i)
        ax1.set_title('Raw data: %04d'%(i))

        ax2 = fig.add_subplot(212, projection='3d')
        ax2.scatter(T[:i+1],Y[:i+1], ac_s[:i+1], label='self', c='blue')
        ax2.scatter(T[:i+1],Y[:i+1], ac_x[:i+1], label='fft', c='orange')
        ax2.set_xlim([0,ac_x.shape[0]])
        ax2.set_ylim([0,ac_x.shape[1]])
        ax2.set_zlim([-1,1])
        ax2.legend()

        ax2.set_title('FFT L:%.3f;\n Self L:%.3f'%(L_fft,L_self))

        plt.savefig("./plots_to_sort/AutocorrelationTest/AT_1D_plots/DD%04d.png"%(i))
        plt.close()
    animator("./plots_to_sort/AutocorrelationTest/AT_1D_plots","AutocorrelationTest_1D")
    

#Test Cases
if 1:
    #0D, 1 number varying through time from -1 to 1, sinusoidal at first then random
    x = np.linspace(0, 2 * np.pi, 50)
    x = np.sin(x)
    x = np.append(x,np.random.uniform(-1, 1, 75))
    plot_autocorrelation_0d(x)
if 1:
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



def Correlation_Plot(path_to_data, number_of_frames):
    #import data
    if not os.path.exists(path_to_output_plots):
        os.makedirs(path_to_output_plots)

    for frame_number in range(0,number_of_frames):
        frame_str = str(frame_number).zfill(4)
        data = yt.load(path_to_data+"DD"+frame_str+"/data"+frame_str)
        all_data = data.all_data()
        print(data.field_list)

        plot = yt.ProjectionPlot(data, coordinate, dataset_to_plot, fontsize = 10)



