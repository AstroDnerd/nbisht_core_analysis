#!/bin/sh
#----------------------------------------------------
# Example SLURM job script to run multiple mpi
# applications within one batch job on TACC's
# Stampede system.
#----------------------------------------------------

#SBATCH -J GPUTrain_convGRUNet_model1_0     # Job name
#SBATCH -o GPUTrain_convGRUNet_model1_0.o%j # Name of stdout output file(%j expands to jobId)
#SBATCH -e GPUTrain_convGRUNet_model1_0.o%j # Name of stderr output file(%j expands to jobId)
#SBATCH -p ai          # Submit to the 'normal' or 'development' queue
#SBATCH -t 24:00:00             # Run time (hh:m) - 48 hours
#SBATCH -N 1                    # Total number of nodes requested (16 cores/node)
#SBATCH -c 48                   #number of CPU cores
#SBATCH --gpus-per-node=4  # Number of GPUs per node
#SBATCH --mail-user=npb22a@fsu.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
#SBATCH -A phy240036
# The next line is required if the user has more than one project

#USAGE 

set echo
source ~/.bashrc

python convGRUNet_train.py