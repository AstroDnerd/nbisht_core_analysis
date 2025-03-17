# standard system modules
import os, sys
os.environ["PATH"] += os.pathsep + "/home/nbisht/myapps/bin/"
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


# standard module for high-quality plots
from PIL import Image
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
mp.rcParams.update(mp.rcParamsDefault)
mp.rcParams['agg.path.chunksize'] = 100000
import matplotlib.gridspec as gridspec
import scienceplots

plt.style.use(['science','grid'])
import sklearn.metrics as skm

# set a seed to ensure reproducibility
seed = 128
rnd  = np.random.RandomState(seed)

MCDSO_CORE = '/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_ML_dataset_XGBoost_MCDS0_core.csv'
MCDSO_NONCORE = '/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_ML_dataset_XGBoost_MCDS0_noncore.csv'
CORESET  = '/data/cb1/nbisht/anvil_scratch/projects/128/B2/datasets/nb101_all_frames.h5'

from matplotlib.collections import LineCollection


def colored_line(x, y, c, ax, **lc_kwargs):
    """
    Plot a line with a color specified along the line by a third value.

    It does this by creating a collection of line segments. Each line segment is
    made up of two straight lines each connecting the current (x, y) point to the
    midpoints of the lines connecting the current point with its two neighbors.
    This creates a smooth line with no gaps between the line segments.

    Parameters
    ----------
    x, y : array-like
        The horizontal and vertical coordinates of the data points.
    c : array-like
        The color values, which should be the same size as x and y.
    ax : Axes
        Axis object on which to plot the colored line.
    **lc_kwargs
        Any additional arguments to pass to matplotlib.collections.LineCollection
        constructor. This should not include the array keyword argument because
        that is set to the color argument. If provided, it will be overridden.

    Returns
    -------
    matplotlib.collections.LineCollection
        The generated line collection representing the colored line.
    """
    if "array" in lc_kwargs:
        warnings.warn('The provided "array" keyword argument will be overridden')

    # Default the capstyle to butt so that the line segments smoothly line up
    default_kwargs = {"capstyle": "butt"}
    default_kwargs.update(lc_kwargs)

    # Compute the midpoints of the line segments. Include the first and last points
    # twice so we don't need any special syntax later to handle them.
    x = np.asarray(x)
    y = np.asarray(y)
    x_midpts = np.hstack((x[0], 0.5 * (x[1:] + x[:-1]), x[-1]))
    y_midpts = np.hstack((y[0], 0.5 * (y[1:] + y[:-1]), y[-1]))

    # Determine the start, middle, and end coordinate pair of each line segment.
    # Use the reshape to add an extra dimension so each pair of points is in its
    # own list. Then concatenate them to create:
    # [
    #   [(x1_start, y1_start), (x1_mid, y1_mid), (x1_end, y1_end)],
    #   [(x2_start, y2_start), (x2_mid, y2_mid), (x2_end, y2_end)],
    #   ...
    # ]
    coord_start = np.column_stack((x_midpts[:-1], y_midpts[:-1]))[:, np.newaxis, :]
    coord_mid = np.column_stack((x, y))[:, np.newaxis, :]
    coord_end = np.column_stack((x_midpts[1:], y_midpts[1:]))[:, np.newaxis, :]
    segments = np.concatenate((coord_start, coord_mid, coord_end), axis=1)

    lc = LineCollection(segments, **default_kwargs)
    lc.set_array(c)  # set the colors of each segment

    return ax.add_collection(lc)


color = np.linspace(50, 121, 71)

def get_df(do_core):
    if do_core == 1:
        df_timeseries = pd.read_csv(MCDSO_CORE)
        model_prefix = 'Core'
    elif do_core == 0:
        df_timeseries = pd.read_csv(MCDSO_NONCORE)
        model_prefix = 'NonCore'
    else:
        df_timeseries_core = pd.read_csv(MCDSO_CORE)
        df_timeseries_noncore = pd.read_csv(MCDSO_CORE)
        df_timeseries = pd.concat([df_timeseries_core, df_timeseries_noncore], axis=0).drop_duplicates()
        model_prefix = 'Combined'
    print(model_prefix)
    return df_timeseries

#get core stuff
with h5py.File(CORESET, 'r') as f:
    tm = f['track_manager']
    all_core_ids = tm['core_ids'][()]
    df_core_pids = tm['particle_ids'][()]

df_core_particles = pd.DataFrame({'Particle_id': df_core_pids, 'Core_id': all_core_ids})


#Same Particle_ids across Initial_Frames
def get_track(df, core_id):
    pid_list = df_core_particles[df_core_particles['Core_id'] == core_id]['Particle_id'].unique()
    print(core_id,len(pid_list))
    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(111)
    if 1000>len(pid_list)>500:
        pid_list = pid_list[::2]
    elif len(pid_list)>1000:
        pid_list = pid_list[::5]
    for particle_id in pid_list:
        track_df = df[df['Particle_id'] == particle_id]
        xpos = track_df['X_f'].values
        ypos = track_df['Y_f'].values
        zpos = track_df['Z_f'].values
        if min(xpos) < 0.01 or max(xpos) > 0.99:
            print('Not This One')
            plt.close()
            return 1
        if min(ypos) < 0.01 or max(ypos) > 0.99:
            plt.close()
            print('Not This One')
            return 1
        lines = colored_line(xpos, ypos, color, ax, cmap="plasma_r", alpha = 0.5)
        ax.scatter(xpos[::10], ypos[::10], c = 'grey', s=5, marker = '.')
    fig.colorbar(lines)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Core : '+str(core_id))
    plt.savefig('./plots_to_sort/tracks_2D/track_2D_core_'+str(core_id)+'.png')
    plt.close()


def plot_features_framewise(df):
    frames = df['Initial_Frame'].unique()
    for frame in frames[0::10]:
        print(frame)
        df_frame = df[df['Initial_Frame'] == frame]
        vx = df_frame['Vx_i'].values
        vy = df_frame['Vy_i'].values
        vz = df_frame['Vz_i'].values
        v = np.sqrt(vx**2 + vy**2 + vz**2)
        rho = df_frame['Density_i'].values
        fig = plt.figure(figsize=(20, 20))
        ax = fig.add_subplot(211)
        ax.set_ylabel('V', fontsize=18)
        ax.hist(v, bins=50, color='black', alpha=0.5)
        ax = fig.add_subplot(212)
        ax.set_ylabel('Density', fontsize=18)
        ax.hist(np.log10(rho), bins=50, color='purple', alpha=0.5)
        
        plt.savefig('./plots_to_sort/Features_framewise/features_frame_'+str(frame)+'.png')
        plt.close()

def plot_features_time(df, threshold = 1000):
    for core_id in np.unique(df_core_particles['Core_id']):
        if len(df_core_particles[df_core_particles['Core_id'] == core_id]['Particle_id'].unique())>threshold:
            fig = plt.figure(figsize=(20, 20))
            ax1 = fig.add_subplot(211)
            ax2 = fig.add_subplot(212, sharex=ax1)
            pid_list = df_core_particles[df_core_particles['Core_id'] == core_id]['Particle_id'].unique()
            print(core_id, len(pid_list))
            v_all = []
            rho_all = []
            for particle_id in pid_list[::10]:
                track_df = df[df['Particle_id'] == particle_id]
                vx = track_df['Vx_i'].values
                vy = track_df['Vy_i'].values
                vz = track_df['Vz_i'].values
                v = np.sqrt(vx**2 + vy**2 + vz**2)
                v_all.append(v)
                rho = np.log10(track_df['Density_i'].values)
                rho_all.append(rho)
                frames = track_df['Initial_Frame'].values
                ax1.plot(frames, v, c = 'grey', alpha=0.5, lw=0.5)
                ax2.plot(frames, rho, c = 'grey', alpha=0.5, lw=0.5)
            v_mean, v_25, v_75 = np.mean(v_all, axis=0), np.percentile(v_all, 25, axis=0), np.percentile(v_all, 75, axis=0)
            rho_mean, rho_25, rho_75 = np.mean(rho_all, axis=0), np.percentile(rho_all, 25, axis=0), np.percentile(rho_all, 75, axis=0)
            ax1.plot(frames, v_mean, c = 'black', lw=2)
            ax1.fill_between(frames, v_25, v_75, color = 'black', alpha = 0.25)
            ax2.plot(frames, rho_mean, c = 'purple', lw=2)
            ax2.fill_between(frames, rho_25, rho_75, color = 'purple', alpha = 0.25)
            ax2.set_xlabel('Frame')
            ax1.set_ylabel('V')
            ax2.set_ylabel('Density')
            ax1.set_title('Core : '+str(core_id))
            plt.savefig('./plots_to_sort/Features_time/Features_time_core_'+str(core_id)+'.png')
            plt.close()


if 0:
    df_timeseries = get_df(1)
    for core_id in np.unique(df_core_particles['Core_id']):
        if len(df_core_particles[df_core_particles['Core_id'] == core_id]['Particle_id'].unique())>100:
            get_track(df_timeseries, core_id)

if 0:
    df_timeseries = get_df(-1)
    plot_features_framewise(df_timeseries)

if 1:
    df_timeseries = get_df(1)
    plot_features_time(df_timeseries, 1000)
