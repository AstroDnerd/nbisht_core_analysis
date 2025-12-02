#!/bin/sh
#----------------------------------------------------
# Example SLURM job script to run multiple mpi
# applications within one batch job on TACC's
# Stampede system.
#----------------------------------------------------

#SBATCH -J Train_Suite_model1_0_TRUN     # Job name
#SBATCH -o Train_Suite_model1_0_TRUN.o%j # Name of stdout output file(%j expands to jobId)
#SBATCH -e Train_Suite_model1_0_TRUN.o%j # Name of stderr output file(%j expands to jobId)
#SBATCH -p wholenode          # Submit to the 'normal' or 'development' queue
#SBATCH -t 48:00:00             # Run time (hh:m) - 48 hours
#SBATCH -N 2                    # MAKE IT 2 or 4 NOW Total number of nodes requested (16 cores/node)
#SBATCH -c 128                   #number of CPU cores
#SBATCH --mail-user=npb22a@fsu.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
#SBATCH -A phy240036
# The next line is required if the user has more than one project

#USAGE 
#Fix the node count and time above.
set echo
source ~/.bashrc


nodes=( $( scontrol show hostnames $SLURM_JOB_NODELIST ) )
nodes_array=($nodes)
head_node=${nodes_array[0]}
head_node_ip=$(srun --nodes=1 --ntasks=1 -w "$head_node" hostname --ip-address)

echo Node IP: $head_node_ip
export LOGLEVEL=INFO

srun torchrun \
--nnodes 2 \
--nproc_per_node 4 \
--rdzv_id $RANDOM \
--rdzv_backend c10d \
--rdzv_endpoint $head_node_ip:29500 \
model1_0_convGRUNet_train.py
