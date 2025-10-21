import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import seaborn as sns
from scipy import fftpack
from scipy import ndimage
from skimage import measure
import warnings
import os
warnings.filterwarnings("ignore")

def make_dir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
def _to_numpy(x):
    """Accept torch tensor or numpy array."""
    try:
        import torch
        if isinstance(x, torch.Tensor):
            return x.detach().cpu().numpy()
    except Exception:
        pass
    return np.asarray(x)


# ---------------------------
# Feature distributions per frame (density and velocity magnitude)
# ---------------------------
def plot_feature_distributions(parent_tensor,
                               frames=(20, 55, 90),
                               channels = {'density':0, 'vx':1,'vy':2,'vz':3},
                               bins=128,
                               fig_file=None):
    """
    Plot histograms of log-density and velocity magnitude for given frames.

    parent_tensor : array_like (T, C, Z, Y, X)
    frames : iterable of integers (frame indices)
    channels : dict mapping 'density','vx','vy','vz' to channel indices
    """
    parent_tensor = _to_numpy(parent_tensor)
    T = parent_tensor.shape[0]
    ch = channels
    fig, axes = plt.subplots(2, len(frames), figsize=(4*len(frames), 6), constrained_layout=True)

    for i, f in enumerate(frames):
        assert 0 <= f < T, "frame index out of range"
        slab = parent_tensor[f]  # (C, Z, Y, X)

        # density histogram (log10)
        density = slab[ch['density']]
        # avoid zeros/negatives - add tiny floor then log10
        eps = 1e-12
        log_density = np.log10(np.clip(density, eps, None).ravel())
        axes[0,i].hist(log_density, bins=bins, density=True, alpha=0.8)
        axes[0,i].set_title(f"Frame {f}: log10(density)")
        axes[0,i].set_xlabel("log10(density)")
        axes[0,i].set_ylabel("probability")

        # velocity magnitude histogram
        vx = slab[ch['vx']]
        vy = slab[ch['vy']]
        vz = slab[ch['vz']]
        vmag = np.sqrt(vx**2 + vy**2 + vz**2).ravel()
        axes[1,i].hist(vmag, bins=bins, density=True, alpha=0.8)
        axes[1,i].set_title(f"Frame {f}: |v|")
        axes[1,i].set_xlabel("|v|")
        axes[1,i].set_ylabel("probability")

    if fig_file:
        fig.savefig(fig_file, dpi=200)
    plt.show()


# ---------------------------
# 3D power spectrum (radial average)
# ---------------------------
def radial_profile_3d(ps3d, center=None, nbins=64):
    """
    compute radially averaged profile of a 3D power spectrum ps3d
    ps3d : (Z,Y,X) real-valued (e.g. |FFT|^2)
    returns: k_bins (bin centers), radial_profile
    """
    nz, ny, nx = ps3d.shape
    z = np.arange(nz) - (center[0] if center is not None else nz//2)
    y = np.arange(ny) - (center[1] if center is not None else ny//2)
    x = np.arange(nx) - (center[2] if center is not None else nx//2)
    zz, yy, xx = np.meshgrid(z,y,x, indexing='ij')
    r = np.sqrt(xx**2 + yy**2 + zz**2).ravel()
    vals = ps3d.ravel()
    rmax = r.max()
    bins = np.linspace(0, rmax, nbins+1)
    inds = np.digitize(r, bins) - 1
    radial = np.zeros(nbins)
    counts = np.zeros(nbins)
    for i in range(nbins):
        mask = inds == i
        if mask.sum()>0:
            radial[i] = vals[mask].mean()
            counts[i] = mask.sum()
        else:
            radial[i] = 0
    k = 0.5*(bins[:-1] + bins[1:])
    return k, radial, counts


def plot_power_spectrum_3d(field3d, dx=1.0, nbins=64, fig_file=None, label=None):
    """
    Compute and plot the radially averaged 3D power spectrum of a scalar field.
    field3d: (Z,Y,X)
    dx: physical voxel size (same units along axes)
    """
    arr = _to_numpy(field3d)
    # subtract mean
    arr = arr - arr.mean()
    # fft, shift
    F = np.fft.fftn(arr)
    ps = np.abs(F)**2
    ps = np.fft.fftshift(ps)
    k_bins, radial, counts = radial_profile_3d(ps, center=(arr.shape[0]//2, arr.shape[1]//2, arr.shape[2]//2), nbins=nbins)
    # physical wavenumber scaling: k ~ bin / (N*dx) but we can show relative k
    plt.loglog(k_bins, radial + 1e-20, label=label)
    plt.xlabel("k (arb units)")
    plt.ylabel("Power")
    if label:
        plt.legend()
    if fig_file:
        plt.savefig(fig_file, dpi=200)
    plt.show()
    return k_bins, radial



# ---------------------------
# Core detection and metrics: simple threshold + labeling
# ---------------------------
'''
def detect_cores_from_density(density3d, threshold, min_voxels=8):
    """
    density3d: (Z,Y,X)
    threshold: value to threshold density (in same units as density3d)
    returns: list of dicts with {'label':int, 'mask':bool-array, 'centroid':(x,y,z) in indices, 'volume':nvoxels, 'peak_density':value}
    """
    bw = density3d > threshold
    labeled = measure.label(bw, connectivity=1)
    props = measure.regionprops(labeled, intensity_image=density3d)
    cores = []
    for p in props:
        if p.area < min_voxels:
            continue
        core = {'label': p.label,
                'mask': (labeled == p.label),
                'centroid': p.centroid,
                'volume': p.area,
                'peak_density': p.max_intensity}
        cores.append(core)
    return cores
'''
def detect_cores_from_density(density3d, threshold, min_voxels=None):
    """
    Purely on the basis of thresholding
    returns: list of dicts with {'label':int, 'mask':bool-array, 'centroid':(x,y,z) in indices, 'volume':nvoxels, 'peak_density':value}
    """
    xthresh, ythresh, zthresh = np.where(density3d > threshold)
    cores = []
    for i in range(len(zthresh)):
        core = {'label': i+1,
                'mask': np.zeros_like(density3d, dtype=bool),
                'centroid': (xthresh[i], ythresh[i], zthresh[i]),
                'volume': 1,
                'peak_density': density3d[xthresh[i], ythresh[i], zthresh[i]]}
        core['mask'][xthresh[i], ythresh[i], zthresh[i]] = True
        cores.append(core)
    return cores


def core_detection_metrics(pred_density, true_density, threshold_pred, threshold_true, min_voxels=8):
    """
    Compute simple core detection metrics: precision, recall, F1 based on overlap between detected core masks.

    pred_density, true_density: (Z,Y,X)
    thresholds: threshold_pred and threshold_true are values to binarize predictions and truth
    """
    pred_cores = detect_cores_from_density(pred_density, threshold_pred, min_voxels=min_voxels)
    true_cores = detect_cores_from_density(true_density, threshold_true, min_voxels=min_voxels)

    # build masks
    pred_mask = np.zeros_like(pred_density, dtype=bool)
    for c in pred_cores:
        pred_mask |= c['mask']
    true_mask = np.zeros_like(true_density, dtype=bool)
    for c in true_cores:
        true_mask |= c['mask']

    tp = np.logical_and(pred_mask, true_mask).sum()
    fp = np.logical_and(pred_mask, ~true_mask).sum()
    fn = np.logical_and(~pred_mask, true_mask).sum()
    precision = tp / (tp + fp + 1e-12)
    recall = tp / (tp + fn + 1e-12)
    f1 = 2 * precision * recall / (precision + recall + 1e-12)

    metrics = {'n_true_cores': len(true_cores),
               'n_pred_cores': len(pred_cores),
               'tp_voxels': int(tp),
               'fp_voxels': int(fp),
               'fn_voxels': int(fn),
               'precision': precision,
               'recall': recall,
               'f1': f1}
    return metrics, pred_cores, true_cores


######################################################
#PAPER PLOTTING UTILITIES
######################################################


import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import scipy.stats as stats

#helper functions 
def velocity_magnitude_field(parent_tensor, frame_idx, ch_vx=1, ch_vy=2, ch_vz=3):
    slab = parent_tensor[frame_idx]
    vx = slab[ch_vx]
    vy = slab[ch_vy]
    vz = slab[ch_vz]
    vmag = np.sqrt(vx**2 + vy**2 + vz**2)
    return vmag

def max_intensity_projection(field3d, axis=0):
    # return MIP along axis
    return np.max(field3d, axis=axis)

def midplane_slice(field3d, axis=0):
    idx = field3d.shape[axis] // 2
    if axis == 0:
        return field3d[idx, :, :]
    elif axis == 1:
        return field3d[:, idx, :]
    else:
        return field3d[:, :, idx]

#log density (clip small)
def safe_log10_density(density, eps=1e-12):
    return np.log10(density)

# ---- 1) overview panel -----------------------------------------------------
def plot_sim_overview(parent_tensor, frame_idx, out_file=None, cmap='viridis'):
    """
    Multi-panel overview: density MIP, velocity MIP, power spectrum (calls plot_power_spectrum_3d),
    and log-density histogram with tail inset.
    """
    density = parent_tensor[frame_idx, 0]  # (Z,Y,X)
    vmag = velocity_magnitude_field(parent_tensor, frame_idx)
    # projections
    dens_mip = max_intensity_projection(np.log10(density), axis=0)  # projection on z axis
    vmag_mip = max_intensity_projection(vmag, axis=0)

    fig = plt.figure(figsize=(12,9))
    gs = gridspec.GridSpec(2,3, width_ratios=[1,1,1], height_ratios=[1,1], hspace=0.2, wspace=0.25)

    ax0 = fig.add_subplot(gs[0,0])
    im0 = ax0.imshow(dens_mip.T, origin='lower')
    ax0.set_title("Density MIP (z-proj)")
    ax0.set_xticks([]); ax0.set_yticks([])
    plt.colorbar(im0, ax=ax0, fraction=0.046, pad=0.01)

    ax1 = fig.add_subplot(gs[0,1])
    im1 = ax1.imshow(vmag_mip.T, origin='lower')
    ax1.set_title("Velocity mag MIP (z-proj)")
    ax1.set_xticks([]); ax1.set_yticks([])
    plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.01)

    ax2 = fig.add_subplot(gs[0,2])
    # call your radial power spectrum function if available:
    try:
        k_bins, radial = plot_power_spectrum_3d(density, dx=1.0, nbins=48, fig_file=None, label=None)
        ax2.loglog(k_bins, radial, label='density P(k)')
        ax2.set_xlabel("k (arb)")
        ax2.set_ylabel("Power")
        ax2.set_title("3D power spectrum (density)")
        ax2.legend()
    except Exception as e:
        ax2.text(0.1,0.5, "Power spectrum failed:\n"+str(e))

    # lower-left: log-density histogram (with tail inset)
    ax3 = fig.add_subplot(gs[1,0:2])
    logrho = safe_log10_density(density.ravel())
    # histogram with KDE / smoothing using gaussian_kde on logrho (for visual only)
    vals, bins, _ = ax3.hist(logrho, bins=120, density=False, alpha=0.7, label='pdf')
    ax3.set_xlabel("log10(density)")
    ax3.set_title("Log-density PDF (global)")

    # inset to show tail (zoom rightmost)
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    axins = inset_axes(ax3, width="30%", height="40%", loc='upper right', borderpad=1)
    right = np.percentile(logrho, 100)
    left = np.percentile(logrho, 99.9)
    axins.hist(logrho, bins=50, density=True)
    axins.set_xlim(3, 5)
    axins.set_ylim(0, 5)
    axins.set_yscale('linear')
    axins.set_title("Tail (zoom)")
    axins.tick_params(axis='both', which='major', labelsize=7)

    # lower-right: show text stats
    ax4 = fig.add_subplot(gs[1,2])
    ax4.axis('off')
    stats_text = [
        f"frame: {frame_idx}",
        f"mean(log10 rho) = {np.mean(logrho):.3f}",
        f"std(log10 rho) = {np.std(logrho):.3f}",
        f"max(log10 rho) = {np.max(logrho):.2f}",
        f"99.99 percentile = {np.percentile(logrho, 99.99):.2f}",
        f"median vmag = {np.median(vmag):.3f}",
        f"rms vmag = {np.sqrt(np.mean(vmag**2)):.3f}"
    ]
    ax4.text(0.05, 0.95, "\n".join(stats_text), va='top', fontsize=11, family='monospace')

    if out_file:
        plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.show()

# ---- 2) density PDF, CCDF, mass fraction -----------------------------------
def plot_density_pdf_ccdf_mass(parent_tensor, frame_idx, nbins=256, out_file=None):
    """
    Plot: (a) PDF of log-density; (b) CCDF voxels (fraction > rho); (c) cumulative mass fraction vs density.
    """
    density = parent_tensor[frame_idx, 0]
    rho = density.ravel()
    rho_pos = np.clip(rho, 1e-12, None)
    logrho = np.log10(rho_pos)

    # PDF (log-binned)
    bins = np.linspace(np.percentile(logrho, 0), np.percentile(logrho, 100), nbins)
    pdf_vals, edges = np.histogram(logrho, bins=bins, density=True)
    centers = 0.5*(edges[:-1] + edges[1:])

    # CCDF (voxel count fraction)
    sorted_rho = np.sort(rho_pos)
    voxels_above = np.arange(sorted_rho.size, 0, -1) / sorted_rho.size
    # mass fraction: sort voxels by density descending and compute cumulative mass
    mass = rho_pos  # assume unit mass per voxel times density; if physical voxel mass scaling required, multiply
    order = np.argsort(mass)[::-1]
    mass_sorted = mass[order]
    cumulative_mass_fraction = np.cumsum(mass_sorted) / (np.sum(mass_sorted)+1e-30)
    density_sorted_desc = mass_sorted  # densities in desc order

    # For plotting mass fraction as function of density, we can sample unique density bins:
    uniq_dens, inds = np.unique(np.round(np.log10(density_sorted_desc), 6), return_index=True)
    sample_mass_frac = cumulative_mass_fraction[inds]

    fig, ax = plt.subplots(1,1, figsize=(6,4))
    ax.plot(centers, pdf_vals, label='pdf (log10 rho)')
    ax.set_xlabel('log10(density)')
    ax.set_ylabel('PDF (arb)')
    ax_twin = ax.twinx()
    # CCDF plotted on twin axis
    # use log-density grid to compute CCDF
    ccdf_vals = [np.mean(rho_pos > 10**(x)) for x in centers]
    ax_twin.semilogy(centers, ccdf_vals, color='C1', label='CCDF (vox > rho)')
    ax_twin.set_ylabel('CCDF (fraction voxels > rho)')

    # add inset showing cumulative mass fraction vs density
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    axins = inset_axes(ax, width="35%", height="35%", loc='upper right')
    # plot cumulative mass fraction vs log10(density) threshold:
    log_dens_thresholds = np.linspace(3, np.max(logrho), 200)
    mass_frac_vs_thresh = [np.sum(mass[rho_pos > 10**th]) / (np.sum(mass)+1e-30) for th in log_dens_thresholds]
    axins.plot(log_dens_thresholds, mass_frac_vs_thresh)
    axins.set_xlabel('log10(density thr)', fontsize=8)
    axins.set_ylabel('mass fraction', fontsize=8)
    axins.tick_params(labelsize=8)

    # legends
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax_twin.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, loc='upper left', fontsize=8)

    plt.title(f"Density PDF + CCDF + mass fraction (frame {frame_idx})")
    if out_file:
        plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.show()

# ---- 3) joint density vs velocity + highlight core voxels -------------------
def plot_density_vs_velocity_with_cores(parent_tensor, frame_idx, detect_threshold_percentile=None, direct_thresh = 1000, min_voxels=8, out_file=None):
    """
    Hexbin/joint plot of log10(density) vs |v|. Overplot core voxels (sampled) as red points.
    """
    density = parent_tensor[frame_idx, 0]
    vmag = velocity_magnitude_field(parent_tensor, frame_idx)
    logrho = safe_log10_density(density)
    rho_flat = density.ravel()
    v_flat = vmag.ravel()

    # detect cores - returns masks and centroids
    if detect_threshold_percentile is None:
        thresh = direct_thresh
    else:
        thresh = np.percentile(density, detect_threshold_percentile)
    cores = detect_cores_from_density(density, threshold=thresh, min_voxels=min_voxels)
    # create mask of voxels belonging to cores (risk memory: use boolean)
    core_mask = np.zeros_like(density, dtype=bool)
    for c in cores:
        core_mask |= c['mask']
    core_mask_flat = core_mask.ravel()

    fig = plt.figure(figsize=(6,5))
    ax = fig.add_subplot(111)
    hb = ax.hexbin(logrho.ravel(), v_flat, gridsize=200, bins='log', cmap='Blues')
    plt.colorbar(hb, ax=ax, label='log(count)')
    ax.set_xlabel('log10(density)')
    ax.set_ylabel('|v|')

    # overlay core voxels (sample a subset if too many)
    core_indices = np.where(core_mask_flat)[0]
    n_overlay = 2000
    if core_indices.size > 0:
        samp = np.random.choice(core_indices, size=min(n_overlay, core_indices.size), replace=False)
        ax.scatter(logrho.ravel()[samp], v_flat[samp], s=8, c='red', alpha=0.6, label='core voxels')
    ax.legend()
    ax.set_title(f"Density vs |v| (frame {frame_idx}) threshold {thresh:.2f}")
    if out_file:
        plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.show()

# ---- 4) core overlay on MIP and orthogonal slices --------------------------------
def plot_cores_overlay(parent_tensor, frame_idx, detect_threshold_percentile=None, direct_thresh = 1000, min_voxels=8, out_file=None):
    """
    Plot: density MIP with core centroids overlaid (colored by peak density). Also show orthogonal slices with core contours.
    """
    density = parent_tensor[frame_idx, 0]
    # detect cores
    if detect_threshold_percentile is None:
        thresh = direct_thresh
    else:
        thresh = np.percentile(density, detect_threshold_percentile)
    cores = detect_cores_from_density(density, threshold=thresh, min_voxels=min_voxels)

    # centroids and peak densities
    centroids = [c['centroid'] for c in cores]
    peakvals = [c['peak_density'] for c in cores]
    if len(centroids) == 0:
        print("No cores found for threshold", thresh)

    fig = plt.figure(figsize=(12,4))
    ax1 = fig.add_subplot(1,3,1)
    mip = max_intensity_projection(np.log10(density), axis=0)
    ax1.imshow(mip.T, origin='lower')
    ax1.set_title('Density MIP (z-proj)')
    # overlay centroids (project z->x,y)
    if len(centroids)>0:
        cent = np.array(centroids)
        # cent[:,1], cent[:,2] are y,x when projecting z
        sc = ax1.scatter(cent[:,2], cent[:,1], c=np.log10(peakvals), s=10, cmap='inferno', edgecolor='white')
        plt.colorbar(sc, ax=ax1, label='log10(peak density)')

    # orthogonal slices at midplane
    midz = density.shape[0]//2
    midy = density.shape[1]//2
    midx = density.shape[2]//2
    ax2 = fig.add_subplot(1,3,2)
    ax2.imshow(density[midz,:,:].T, origin='lower')
    ax2.set_title(f'z-slice (z={midz})')
    # overlay core contours (where mask intersects slice)
    for c in cores:
        mask = c['mask']
        if mask[midz].any():
            # get contour of that mask in slice (simple outline via bounding box)
            ys, xs = np.where(mask[midz])
            ax2.scatter(xs, ys, s=8, c='red', alpha=0.4)

    ax3 = fig.add_subplot(1,3,3)
    ax3.imshow(density[:,midy,:].T, origin='lower')
    ax3.set_title(f'y-slice (y={midy})')
    for c in cores:
        # mask transpose coords
        mask = c['mask']
        if mask[:,midy,:].any():
            zs, xs = np.where(mask[:,midy,:])
            ax3.scatter(xs, zs, s=8, c='red', alpha=0.4)

    if out_file:
        plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.show()

# ---- 5) 3D scatter of core centroids (matplotlib) ----------------------------
def plot_core_centroids_3d(parent_tensor, frame_idx, detect_threshold_percentile=None, direct_thresh = 1000, min_voxels=8, out_file=None):
    density = parent_tensor[frame_idx, 0]
    if detect_threshold_percentile is None:
        thresh = direct_thresh
    else:
        thresh = np.percentile(density, detect_threshold_percentile)
    cores = detect_cores_from_density(density, threshold=thresh, min_voxels=min_voxels)
    centroids = np.array([c['centroid'] for c in cores])
    peakvals = np.array([c['peak_density'] for c in cores])
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111, projection='3d')
    if centroids.size == 0:
        ax.text(0.5,0.5,0.5, "No cores", transform=ax.transAxes)
    else:
        sc = ax.scatter(centroids[:,0], centroids[:,1], centroids[:,2],
                        c=np.log10(peakvals), s=40, cmap='plasma', depthshade=True)
        plt.colorbar(sc, ax=ax, label='log10(peak density)')
        ax.set_xlabel('z'); ax.set_ylabel('y'); ax.set_zlabel('x')
        ax.set_title(f'Core centroids (frame {frame_idx}), {len(cores)} cores')
    if out_file:
        plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.show()

def getcoresraw(simname, basepath, ddnumber, coreshdf5_name, c_min = None, c_max = None, step = None):
    import yt
    from yt.data_objects.level_sets.api import Clump,find_clumps,return_bottom_clumps
    YT_x = ('gas','x')
    YT_y = ('gas','y')
    YT_z = ('gas','z')
    simdir = basepath + simname + "/"
    if os.path.exists(coreshdf5_name):
        print("File exists, returning File data.", coreshdf5_name)
        import h5py
        fptr=h5py.File(coreshdf5_name,'r')
        peak_list = fptr['peaks'][:]
        dens_values = fptr['peak_density'][:]
        fptr.close()
        return peak_list, dens_values
    
    ds = yt.load("%s/TT%04d/time%04d"%(simdir,ddnumber,ddnumber))

    #get the clumps
    ad = ds.all_data()
    master_clump = Clump(ad,('gas','density'))
    master_clump.add_validator("min_cells", 8)
    if c_min is None:
        c_min = 10 #ad["gas", "density"].min()
    #c_max = 534069645. # ad["gas", "density"].max()
    if c_max is None:
        c_max = ad["gas", 'density'].max()
    if step is None:
        step = 100
    find_clumps(master_clump, c_min, c_max, step)
    #get the peaks from the clump.  Save in h5name
    leaf_clumps = return_bottom_clumps(master_clump)

    peak_list=[]
    den_max=[]
    for i in range(len(leaf_clumps)):
        den_max.append(leaf_clumps[i][('gas','density')].max())
        max_loc = np.where(leaf_clumps[i]['gas','density']==den_max[i])
        a = leaf_clumps[i][YT_x][max_loc][0]
        b = leaf_clumps[i][YT_y][max_loc][0]
        c = leaf_clumps[i][YT_z][max_loc][0]

        this_peak = ds.arr([a,b,c],'code_length')
        peak_list.append(this_peak)
    dens_values = ds.find_field_values_at_points(('gas','density'), peak_list)
    dens_values = np.array(dens_values, dtype='float64')
    import h5py
    fptr=h5py.File(coreshdf5_name,'w')
    fptr.create_dataset('peaks',data=np.array(peak_list))
    fptr.create_dataset('peak_density',data=dens_values)
    fptr.close()
    print("made peaks in %s"%coreshdf5_name)

    return peak_list, dens_values

simdata = []
if __name__ == "__main__":
    #parent_tensor shape (100,4,128,128,128)
    basepath = "/anvil/scratch/x-nbisht1/projects/512/NonsinkSimSuite/"
    path_to_output_plots = "/home/x-nbisht1/projects/results/512/NonsinkSimSuite/CubeStats/"
    simarray = ['d03_Ms1.0_Ma0.0_512', 'd07_Ms2.0_Ma0.0_512', 'd11_Ms3.0_Ma0.0_512', 
            'd15_Ms4.0_Ma0.0_512', 'd19_Ms5.0_Ma0.0_512', 'd23_Ms6.0_Ma0.0_512', 
            'd27_Ms7.0_Ma0.0_512', 'd31_Ms8.0_Ma0.0_512']
    ddnumber = ['DD0030', 'DD0050', 'DD0070', 'DD0090']
    make_dir(path_to_output_plots)
    for sim in simarray[0:1]:
        for dd in ddnumber[0:1]:
            datasetname = sim[0:3]+'_'+dd #'d03_DD0030'
            simname = sim+'_'+dd #'d03_Ms1.0_Ma0.0_512_DD0030'
            TSCube = []
            for fieldname in ['density', 'x-velocity', 'y-velocity', 'z-velocity']:
                df_name = '/anvil/scratch/x-nbisht1/projects/512/NonsinkSimSuite/timecubes/'+datasetname+'_TimeseriesCubes_'+fieldname+'.npy'
                infile = open(df_name, 'rb')
                TSCube_field = np.load(infile)
                TSCube.append(TSCube_field)
                time_array = np.load(infile)
                DELTAT = np.mean(time_array[1:]-time_array[:-1])
                infile.close()
            TSCube = np.stack(TSCube, axis=1) #(100,4,128,128,128)
            #plot feature histograms at frames 20, 55, 90
            if 0:
                plot_feature_distributions(TSCube, frames=(20,55,90), fig_file=path_to_output_plots+"/feat_hist_"+datasetname+".png")
                #compute 3D power spectrum of one snapshot density
                snapshot = TSCube[55,0]  # density at frame 55
                plot_power_spectrum_3d(snapshot, dx=1.0, nbins=64, fig_file=path_to_output_plots+"/power_spec55_"+datasetname+".png")
                thresh = np.percentile(snapshot, 95) #95th percentile threshold
                pixel_cores = detect_cores_from_density(snapshot, threshold = thresh, min_voxels=8)
                print(f"{datasetname}: detected {len(pixel_cores)} cores at 95th percentile threshold")
                print(pixel_cores)
            if 1:
                frame_idx = 99
                plot_sim_overview(TSCube, frame_idx, out_file=path_to_output_plots+"/sim_overview_"+datasetname+".png", cmap='viridis')
                plot_density_pdf_ccdf_mass(TSCube, frame_idx, nbins=256, out_file=path_to_output_plots+"/density_pdf_ccdf_mass_"+datasetname+".png")
                #detect_threshold_percentile=99.9
                #detect_threshold_percentile_str = str(detect_threshold_percentile).replace('.','p')
                detect_threshold_percentile = None
                thresh=2500
                detect_threshold_percentile_str = '10pow3'
                min_voxels=1
                plot_density_vs_velocity_with_cores(TSCube, frame_idx, detect_threshold_percentile=detect_threshold_percentile, direct_thresh = thresh, min_voxels=min_voxels,
                                                     out_file=path_to_output_plots+"/density_vs_velocity_with_cores_"+datasetname+
                                                     "_thresh_"+detect_threshold_percentile_str+"_vox_"+str(min_voxels)+".png")
                plot_cores_overlay(TSCube, frame_idx, detect_threshold_percentile=detect_threshold_percentile, direct_thresh = thresh, min_voxels=min_voxels,
                                                     out_file=path_to_output_plots+"/cores_overlay_"+datasetname+
                                                     "_thresh_"+detect_threshold_percentile_str+"_vox_"+str(min_voxels)+".png")
                plot_core_centroids_3d(TSCube, frame_idx, detect_threshold_percentile=detect_threshold_percentile, direct_thresh = thresh, min_voxels=min_voxels,
                                                     out_file=path_to_output_plots+"/core_centroids_3d_"+datasetname+
                                                     "_thresh_"+detect_threshold_percentile_str+"_vox_"+str(min_voxels)+".png")
                coreshdf5_name = '/anvil/scratch/x-nbisht1/projects/512/NonsinkSimSuite/dataset/'+datasetname+'_cores.hdf5'
                peak_list, density_vals = getcoresraw(simname, basepath, frame_idx, coreshdf5_name, c_min = None, c_max = None, step = None)
                simdata.append(peak_list)
                dens_mask = density_vals > thresh
                cores_raw = np.array(peak_list)[dens_mask]
                cores_density = density_vals[dens_mask]
                print(f"{datasetname}: detected {len(cores_raw)} cores from yt clump finder at threshold {thresh}")
                #compare core locations
                cores = detect_cores_from_density(TSCube[frame_idx, 0], threshold=thresh, min_voxels=min_voxels)
                cores_dataset = []
                cores_dataset_peakdens = []
                for c in cores:
                    xc,yc,zc = c['centroid']
                    cores_dataset.append([(xc+0.5)/128, (yc+0.5)/128, (zc+0.5)/128])  #normalize to [0,1]
                    cores_dataset_peakdens.append(c['peak_density'])
                cores_dataset = np.array(cores_dataset)
                cores_dataset_peakdens = np.array(cores_dataset_peakdens)
                print("Pixel-based detected cores:", len(cores_dataset))
                #now compare with cores_raw using cKDTree
                if len(cores_raw)>0 and len(cores_dataset)>0:
                    matched_dataset_indices = []
                    from scipy.spatial import cKDTree
                    tree = cKDTree(cores_dataset)
                    dists, indices = tree.query(cores_raw, k=1)
                    match_radius = 1/16  #min distance for correct prediction
                    n_matches = np.sum(dists < match_radius)
                    print(f"Number of matching cores within {match_radius}: {n_matches} out of {len(cores_dataset)} (unit) and {len(cores_raw)} (yt). ")
                    #print matched core details
                    for i, dist in enumerate(dists):
                        if dist < match_radius:
                            idx_pixel = indices[i]
                            print(f"Matched core {i} (yt) at {cores_raw[i]} with core {idx_pixel} (unit) at {cores_dataset[idx_pixel]}, distance = {dist:.5f}, "
                                  f"yt peak density = {cores_density[i]:.2e}, unit peak density = {cores_dataset_peakdens[idx_pixel]:.2e}")
                            matched_dataset_indices.append(idx_pixel)
                        else:
                            print(f"No match for core {i} (yt) at position {cores_raw[i]}, closest distance = {dist:.5f} and density = {cores_density[i]:.2e}")
                    #also print unmatched cores from dataset
                    matched_dataset_indices = set(matched_dataset_indices)
                    for j in range(len(cores_dataset)):
                        matched_cores = cores_dataset[list(matched_dataset_indices)]
                        if j not in matched_dataset_indices:
                            #check if this unmatched core is adjacent to any matched core
                            xi,yi,zi = cores_dataset[j]
                            if any(np.linalg.norm(matched_cores - np.array([xi,yi,zi]), axis=1) < (1/64)):
                                print(f"Core {j} (unit) at position {cores_dataset[j]} is adjacent to a matched core, peak density = {cores_dataset_peakdens[j]:.2e}")
                                matched_dataset_indices.add(j)
                                continue
                            print(f"Unmatched core {j} (unit) at position {cores_dataset[j]}, peak density = {cores_dataset_peakdens[j]:.2e}")
                else:
                    print("No cores to compare between methods.")

    pass
