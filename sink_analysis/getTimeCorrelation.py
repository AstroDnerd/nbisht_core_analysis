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


#define paths sinkoff
path_to_data_128 = "/data/cb1/nbisht/anvil_scratch/projects/128/B2/"
path_to_output_plots = "./plots_to_sort/TimeCorrelation/"


def ac_fft(x):
    N = len(x)
    x_fft = fft.fft(x,n=2*N)
    x_fft_conj = np.conj(x_fft)
    autocorr = fft.ifft(x_fft*x_fft_conj)
    autocorr = autocorr[:N]
    return np.real(autocorr/autocorr[0])

def ac_self(x):
    N = len(x)
    ac_self = []
    denominator = np.sum(np.multiply(x,x))/len(x)
    for i in range(N):
        pdt_arr = np.multiply(x[i:],x[:N-i])
        ac_self.append((np.sum(pdt_arr)/len(pdt_arr))/denominator)
    ac_self = np.array(ac_self)
    return ac_self

x = np.linspace(0, 2 * np.pi, 100)
x = np.sin(x)
x = np.append(x,np.random.rand(25))
ac_x = ac_fft(x)
ac_s = ac_self(x)
print(ac_x)
print(ac_s)
print(np.sum(ac_x))
print(np.sum(ac_s)

plt.plot(ac_x, label='fft')
plt.plot(ac_s, label='self')
plt.legend()
plt.savefig("./plots_to_sort/AutocorrelationTest.png")


def make_dir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

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



