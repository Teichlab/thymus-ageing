#!/bin/bash
#SBATCH --job-name=scvi_modality
#SBATCH --account=TEICHLAB-SL2-GPU
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=20:00:00
#SBATCH --output=logs/scvi_%j.out
#SBATCH --error=logs/scvi_%j.err

set -eo pipefail

eval "$(conda shell.bash hook)"
conda activate scvi_env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

python /home/nd590/rfs/rfs-teichlab-nfs-iCNyzSAaucw/nd590/projects/aging/scripts_spatial/03_run_scvi_aging.py
