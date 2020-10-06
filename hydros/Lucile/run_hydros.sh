#!/bin/bash -l
#SBATCH --time=02:00:00
#SBATCH --job-name=run_hydros
#SBATCH --partition=short
#SBATCH -n 1
##SBATCH --mem=25G

srun python run_hydros.py
