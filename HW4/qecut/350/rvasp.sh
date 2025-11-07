#!/bin/bash
# --- SLURM Resource Settings ---
#SBATCH --job-name=VASP_run    # A name for your job
#SBATCH --output=slurm.out   # File to save standard output
#SBATCH --error=slurm.err    # File to save errors
#SBATCH --nodes=1            # Request 1 node
#SBATCH --ntasks=1           # Request 1 CPU core
#SBATCH --partition=parallel # Partition you found earlier
#SBATCH --time=01:00:00      # 1 hour of run time

# --- Email Notifications ---
#SBATCH --mail-user=mkarumu1@jh.edu   # Your email address
#SBATCH --mail-type=ALL               # Notify on BEGIN, END, and FAIL

# --- Your Commands ---
# 1. Set up the environment (This line is essential)
# source ~/scr4-en510/setup_Environment.sh

# 2. Run your VASP calculation
mpirun -n 1 vasp_std