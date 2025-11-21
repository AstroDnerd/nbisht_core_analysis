import torch
import torch.profiler
import time
import cProfile
import pstats
import io
from contextlib import contextmanager
from collections import defaultdict
import json
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

class TrainingProfiler:
    """Comprehensive profiler for training pipeline"""
    
    def __init__(self, output_dir="profiling_results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.timings = defaultdict(list)
        self.memory_stats = []
        
    @contextmanager
    def timer(self, name):
        """Context manager for timing code sections"""
        start = time.perf_counter()
        try:
            yield
        finally:
            elapsed = time.perf_counter() - start
            self.timings[name].append(elapsed)
            print(f"Timer {name}: {elapsed:.4f}s")
    
    def record_memory(self, label):
        """Record GPU memory usage"""
        if torch.cuda.is_available():
            allocated = torch.cuda.memory_allocated() / 1024**3  # GB
            reserved = torch.cuda.memory_reserved() / 1024**3
            self.memory_stats.append({
                'label': label,
                'allocated_gb': allocated,
                'reserved_gb': reserved,
                'time': time.perf_counter()
            })
            print(f"Record_Memory {label}: {allocated:.2f}GB allocated, {reserved:.2f}GB reserved")
    
    def summarize(self):
        """Print summary of all timings"""
        print("\n" + "="*80)
        print("PROFILING SUMMARY")
        print("="*80)
        
        total_time = 0
        timing_data = []
        
        for name, times in sorted(self.timings.items()):
            avg = np.mean(times)
            std = np.std(times)
            total = np.sum(times)
            count = len(times)
            total_time += total
            
            timing_data.append({
                'name': name,
                'avg': avg,
                'std': std,
                'total': total,
                'count': count
            })
            
            print(f"\n{name}:")
            print(f"  Count: {count}")
            print(f"  Average: {avg:.4f}s")
            print(f"  Std Dev: {std:.4f}s")
            print(f"  Total: {total:.4f}s")
        
        print(f"\n{'='*80}")
        print(f"TOTAL TIME: {total_time:.4f}s")
        print(f"{'='*80}\n")
        
        return timing_data
    
    def plot_results(self):
        """Create visualization of profiling results"""
        if not self.timings:
            print("No timing data to plot")
            return
        
        # Plot 1: Time breakdown pie chart
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # Aggregate timings by category
        categories = {
            'Data Loading': ['load_batch', 'data_transfer', 'batch_unpack'],
            'Forward Pass': ['forward_pass', 'model_forward'],
            'Loss Computation': ['loss_computation', 'compute_losses'],
            'Backward Pass': ['backward_pass', 'gradient_step'],
            'Other': []
        }
        
        category_times = defaultdict(float)
        for name, times in self.timings.items():
            total = sum(times)
            categorized = False
            for cat, keywords in categories.items():
                if any(kw in name.lower() for kw in keywords):
                    category_times[cat] += total
                    categorized = True
                    break
            if not categorized:
                category_times['Other'] += total
        
        # Pie chart
        axes[0, 0].pie(category_times.values(), labels=category_times.keys(), 
                       autopct='%1.1f%%', startangle=90)
        axes[0, 0].set_title('Time Distribution by Category')
        
        # Bar chart of individual operations
        names = list(self.timings.keys())[:15]  # Top 15
        totals = [sum(self.timings[name]) for name in names]
        axes[0, 1].barh(names, totals)
        axes[0, 1].set_xlabel('Total Time (s)')
        axes[0, 1].set_title('Top 15 Operations by Total Time')
        axes[0, 1].invert_yaxis()
        
        # Memory usage over time
        if self.memory_stats:
            mem_times = [s['time'] - self.memory_stats[0]['time'] for s in self.memory_stats]
            mem_alloc = [s['allocated_gb'] for s in self.memory_stats]
            mem_reserved = [s['reserved_gb'] for s in self.memory_stats]
            
            axes[1, 0].plot(mem_times, mem_alloc, 'b-', label='Allocated', linewidth=2)
            axes[1, 0].plot(mem_times, mem_reserved, 'r--', label='Reserved', linewidth=2)
            axes[1, 0].set_xlabel('Time (s)')
            axes[1, 0].set_ylabel('Memory (GB)')
            axes[1, 0].set_title('GPU Memory Usage Over Time')
            axes[1, 0].legend()
            axes[1, 0].grid(True, alpha=0.3)
        
        # Batch processing time histogram
        if 'batch_processing' in self.timings:
            axes[1, 1].hist(self.timings['batch_processing'], bins=30, edgecolor='black')
            axes[1, 1].set_xlabel('Time (s)')
            axes[1, 1].set_ylabel('Frequency')
            axes[1, 1].set_title('Batch Processing Time Distribution')
            axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        output_path = self.output_dir / 'profiling_summary.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"\n Plots saved to: {output_path}")
        plt.close()
    
    def save_json(self):
        """Save profiling data to JSON"""
        data = {
            'timings': {k: {'values': v, 'mean': np.mean(v), 'sum': np.sum(v)} 
                       for k, v in self.timings.items()},
            'memory_stats': self.memory_stats
        }
        output_path = self.output_dir / 'profiling_data.json'
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"JSON data saved to: {output_path}")


def profile_data_loading(hdf5_path, num_samples, test_percentage, seed, profiler):
    """Profile HDF5 data loading"""
    print("\n" + "="*80)
    print("PROFILING DATA LOADING")
    print("="*80)
    
    with profiler.timer("total_data_loading"):
        with profiler.timer("hdf5_open"):
            import h5py
            f = h5py.File(hdf5_path, 'r')
        
        with profiler.timer("dataset_creation"):
            # Import your data loading function
            from all_modules import load_hdf5_data_for_training
            train_dataset, test_dataset, train_indices, test_indices = \
                load_hdf5_data_for_training(hdf5_path, num_samples=num_samples, 
                                           test_percentage=test_percentage, seed=seed)
    
    profiler.record_memory("after_data_loading")
    return train_dataset, test_dataset, train_indices, test_indices


def profile_model_initialization(modelclass, argsUNET, argsGRU, device, profiler):
    """Profile model initialization"""
    print("\n" + "="*80)
    print("PROFILING MODEL INITIALIZATION")
    print("="*80)
    
    with profiler.timer("model_creation"):
        model = modelclass(argsUNET, argsGRU)
    
    with profiler.timer("model_to_device"):
        model = model.to(device)
    
    profiler.record_memory("after_model_init")
    
    # Count parameters
    total_params = sum(p.numel() for p in model.parameters())
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"Total parameters: {total_params:,}")
    print(f"Trainable parameters: {trainable_params:,}")
    
    return model


def profile_single_batch(model, batch_data, config, device, profiler, use_amp=True):
    """Profile a single batch in detail"""
    
    with profiler.timer("batch_unpack"):
        if len(batch_data) == 3:
            x, y, labels = batch_data
        else:
            x, y = batch_data
    
    with profiler.timer("data_transfer"):
        x = x.to(device, non_blocking=True)
        y = y.to(device, non_blocking=True)
    
    profiler.record_memory("after_data_transfer")
    
    with profiler.timer("forward_pass"):
        if use_amp and device.type == 'cuda':
            with torch.amp.autocast(device_type='cuda'):
                pred = model(x)
            pred_f = pred.float()
            y_f = y.float()
        else:
            pred = model(x)
            pred_f = pred
            y_f = y
    
    profiler.record_memory("after_forward")
    
    with profiler.timer("loss_computation"):
        from all_modules import return_total_loss_multichannel
        l1, hist_l, spectral_l, hd_l, mass_l = return_total_loss_multichannel(
            y_f, pred_f, config
        )
        losses = torch.stack([l1, hist_l, mass_l, spectral_l, hd_l])
        total_loss = (torch.tensor(config.get('static_w', [1,1,1,1,1]), 
                                  device=device) * losses).sum()
    
    profiler.record_memory("after_loss")
    
    with profiler.timer("backward_pass"):
        if use_amp and device.type == 'cuda':
            scaler = torch.amp.GradScaler()
            scaler.scale(total_loss).backward()
        else:
            total_loss.backward()
    
    profiler.record_memory("after_backward")
    
    return total_loss.item()


def run_pytorch_profiler(model, loader, device, output_dir, num_batches=5):
    """Run PyTorch's built-in profiler for detailed analysis"""
    print("\n" + "="*80)
    print("RUNNING PYTORCH PROFILER")
    print("="*80)
    
    activities = [torch.profiler.ProfilerActivity.CPU]
    if device.type == 'cuda':
        activities.append(torch.profiler.ProfilerActivity.CUDA)
    
    with torch.profiler.profile(
        activities=activities,
        record_shapes=True,
        profile_memory=True,
        with_stack=True,
        with_flops=True,
        with_modules=True,
        on_trace_ready=torch.profiler.tensorboard_trace_handler(str(output_dir / 'torch_profiler'))
    ) as prof:
        model.train()
        for batch_idx, batch_data in enumerate(loader):
            if batch_idx >= num_batches:
                break
            
            x, y = batch_data[0], batch_data[1]
            x, y = x.to(device), y.to(device)
            
            pred = model(x)
            loss = torch.nn.functional.mse_loss(pred, y)
            loss.backward()
            
            prof.step()
    
    # Print summary
    print("\n" + "-"*80)
    print("Top 10 operations by CPU time:")
    print("-"*80)
    print(prof.key_averages().table(sort_by="cpu_time_total", row_limit=10))
    
    if device.type == 'cuda':
        print("\n" + "-"*80)
        print("Top 10 operations by CUDA time:")
        print("-"*80)
        print(prof.key_averages().table(sort_by="cuda_time_total", row_limit=10))
    
    # Save detailed results
    #prof.export_chrome_trace(str(output_dir / "torch_trace.json"))
    print(f"\n Chrome trace saved to: {output_dir / 'torch_trace.json'}")
    print("   View at: chrome://tracing")


def main_profiling_pipeline(master_script_params):
    """Main profiling function that wraps your training pipeline"""
    
    profiler = TrainingProfiler()
    
    # Extract parameters
    DATASETNAME = master_script_params['DATASETNAME']
    number_of_samples = master_script_params['number_of_samples']
    TEST_PERCENTAGE = master_script_params['TEST_PERCENTAGE']
    seed = master_script_params['seed']
    DEVICE = master_script_params['DEVICE']
    argsUNET = master_script_params['argsUNET']
    argsGRU = master_script_params['argsGRU']
    modelclass = master_script_params['modelclass']
    
    # Profile data loading
    train_dataset, test_dataset, train_indices, test_indices = profile_data_loading(
        DATASETNAME, number_of_samples, TEST_PERCENTAGE, seed, profiler
    )
    
    # Profile model initialization
    model = profile_model_initialization(modelclass, argsUNET, argsGRU, DEVICE, profiler)
    
    # Create dataloader
    with profiler.timer("dataloader_creation"):
        from torch.utils.data import DataLoader
        loader = DataLoader(
            train_dataset, 
            batch_size=argsUNET.batch_size,
            shuffle=True,
            num_workers=argsUNET.num_workers,
            pin_memory=(DEVICE.type == 'cuda'),
            persistent_workers=(argsUNET.num_workers > 0)
        )
    
    # Profile first few batches in detail
    print("\n" + "="*80)
    print("PROFILING INDIVIDUAL BATCHES")
    print("="*80)
    
    # Setup config for loss computation
    config = {
        'channel_scales': [1.0, 1.0, 1.0, 1.0],
        'channel_weights': [1.0, 1.0, 1.0, 1.0],
        'density_channel': 0,
        'velocity_channels': [1,2,3],
        'density_is_log': True,
        'hist_bins': 32,
        'hist_vmin': -5.0,
        'hist_vmax': 5.0,
        'hist_sigma': None,
        'hist_sample_frac': 0.2,
        'do_spectral': True,
        'spec_sample_frac': 1.0,
        'hd_q': 0.999,
        'hd_smooth': 1e-2,
        'static_w': [1.0, 1.0, 1.0, 1.0, 1.0]
    }
    
    model.train()
    optimizer = torch.optim.AdamW(model.parameters(), lr=argsUNET.lr)
    
    for batch_idx, batch_data in enumerate(loader):
        if batch_idx >= 3:  # Profile first 3 batches
            break
        
        print(f"\n--- Batch {batch_idx + 1} ---")
        with profiler.timer(f"batch_processing"):
            optimizer.zero_grad()
            loss = profile_single_batch(
                model, batch_data, config, DEVICE, profiler,
                use_amp=(DEVICE.type == 'cuda')
            )
            with profiler.timer("optimizer_step"):
                optimizer.step()
        
        if DEVICE.type == 'cuda':
            torch.cuda.synchronize()
    
    # Run PyTorch's built-in profiler
    run_pytorch_profiler(model, loader, DEVICE, profiler.output_dir, num_batches=5)
    
    # Generate summary and visualizations
    profiler.summarize()
    profiler.plot_results()
    profiler.save_json()
    
    print("\n" + "="*80)
    print("PROFILING COMPLETE!")
    print("="*80)
    print(f"\nResults saved to: {profiler.output_dir}")
    print("\nNext steps:")
    print("1. Check profiling_summary.png for visual breakdown")
    print("2. Review profiling_data.json for detailed timing data")
    print("3. Open torch_trace.json in chrome://tracing for detailed GPU analysis")
    print("4. Check models/logs/torch_profiler with TensorBoard")
    
    return profiler

