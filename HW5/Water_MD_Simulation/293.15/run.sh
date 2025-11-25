#!/bin/bash
#SBATCH --job-name=293.15_water_sim
#SBATCH --time=1-00:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mkarumu1@jh.edu

# Load environment
source ~/scr4-en510/setup_environment.sh

# Configure Intel MPI to work with SLURM
export I_MPI_PMI_LIBRARY=/cm/shared/apps/slurm/current/lib64/libpmi.so
export I_MPI_HYDRA_BOOTSTRAP=slurm

echo "Starting Water simulation..."
srun -n 4 lmp -in water.in
echo "Simulation finished."