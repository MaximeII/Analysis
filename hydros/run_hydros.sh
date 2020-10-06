#!/bin/bash -l
#SBATCH --time=00:30:00
#SBATCH --job-name=run_hydros
#SBATCH --partition=short
#SBATCH -n 1
##SBATCH --mem=25G

python run_hydros.py BCQ
python run_hydros.py BGA
python run_hydros.py BCG

echo Job $SLURM_ARRAY_TASK_ID complete.

