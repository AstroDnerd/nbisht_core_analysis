#!/bin/sh
#----------------------------------------------------
# Example SLURM job script to run multiple mpi
# applications within one batch job on TACC's
# Stampede system.
#----------------------------------------------------

#SBATCH -J Train_UNet_2D_static_attention     # Job name
#SBATCH -o Train_UNet_2D_static_attention.o%j # Name of stdout output file(%j expands to jobId)
#SBATCH -e Train_UNet_2D_static_attention.o%j # Name of stderr output file(%j expands to jobId)
#SBATCH -p highmem          # Submit to the 'normal' or 'development' queue
#SBATCH -t 24:00:00             # Run time (hh:m) - 1.5 hours
#SBATCH -N 1                    # Total number of nodes requested (16 cores/node)
#SBATCH -c 12                   #number of CPUs
#SBATCH --mail-user=npb22a@fsu.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
#SBATCH -A phy240036
# The next line is required if the user has more than one project

#USAGE 
#Fix the node count and time above.


set echo
source ~/.bashrc

python Conditional_Unet_2D_Model_attention.py