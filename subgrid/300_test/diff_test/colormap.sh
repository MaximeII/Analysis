#!/bin/bash
#SBATCH --job-name=colormap
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH -n 1                  ## number of tasks
#SBATCH --mem=25G
#SBATCH --exclude=vorna-436

t=8                   #threads per process

module load GCCcore/8.3.0
module load OpenMPI
module load Eigen

export PTNONINTERACTIVE=1

umask 007
ulimit -c unlimited

#--------------------------------------------------------------------
#---------------------DO NOT TOUCH-----------------------------------
nodes=$SLURM_NNODES
#Vorna has 2 x 8 cores
cores_per_node=16
# Hyperthreading
ht=1
#Change PBS parameters above + the ones here
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')
export OMP_NUM_THREADS=$t
export OMPI_MCA_btl_openib_allow_ib=1
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

srun python colormap_and_VDFs_Box.py 900 1 330
