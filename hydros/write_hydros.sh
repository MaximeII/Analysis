#!/bin/bash -l
#SBATCH --time=00:15:00
#SBATCH --job-name=write_hydros
#SBATCH --partition=short
#SBATCH -n 1
##SBATCH --mem=25G

python write_hydros.py BCQ
python write_hydros.py BGA
python write_hydros.py BCG

python write_hydros_mirror.py BCQ
python write_hydros_mirror.py BGA
python write_hydros_mirror.py BCG

echo Job $SLURM_ARRAY_TASK_ID complete.
