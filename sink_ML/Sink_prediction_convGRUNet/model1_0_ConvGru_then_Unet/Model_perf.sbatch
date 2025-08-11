#!/bin/sh
#----------------------------------------------------
# Example SLURM job script to run multiple mpi
# applications within one batch job on TACC's
# Stampede system.
#----------------------------------------------------

#SBATCH -J convGRUNet_model1_0_histperformance     # Job name
#SBATCH -o convGRUNet_model1_0_histperformance.o%j # Name of stdout output file(%j expands to jobId)
#SBATCH -e convGRUNet_model1_0_histperformance.o%j # Name of stderr output file(%j expands to jobId)
#SBATCH -p highmem          # Submit to the 'normal' or 'development' queue
#SBATCH -t 48:00:00             # Run time (hh:m) - 48 hours
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

python convGRUNet_Model_Metric.py