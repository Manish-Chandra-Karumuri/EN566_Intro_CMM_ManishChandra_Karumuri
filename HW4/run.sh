#!/bin/bash

# --- SLURM Directives ---
# These directives request resources from the cluster. Adjust them based on the complexity
# of your ising.py script (e.g., if it runs longer or needs more memory).

# Set the job name
#SBATCH --job-name=Ising_Part2
# Specify the output file for standard output (will be ising_job_12345.out)
# The %j variable is replaced with the job ID
#SBATCH --output=ising_job_%j.out
# Specify the error file for standard error
#SBATCH --error=ising_job_%j.err
# Set the maximum running time (e.g., 30 minutes)
#SBATCH --time=00:30:00
#SBATCH --partition=shared
# Request 1 task (process)
#SBATCH --ntasks=1
# Request 1 CPU core for the task
#SBATCH --cpus-per-task=1
# Request 4 Gigabytes of memory for the job
#SBATCH --mem=2G
# Email address for notification
#SBATCH --mail-user=mkarumu1@jh.edu
# Notify me when job finishes (END)
#SBATCH --mail-type=ALL

# --- Execution Environment Setup ---



module load python
module load matplotlib

# 2. Run the main command
echo "Executing: python ising.py --part=2"
python ising.py --part=2

echo "Job finished on $(date)"