#!/bin/bash

#SBATCH --nodes 1
#SBATCH --job-name=hex_lattice
#SBATCH --output=/data/slurm.out
#SBATCH --error=/data/slurm.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
srun ./compile.sh --oversubscribe
