#!/bin/bash -l
#SBATCH --job-name=diffusion
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH -n 1
##SBATCH --mem=25G
#SBATCH --exclude=vorna-436

srun python write_hydros.py

