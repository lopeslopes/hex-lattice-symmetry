#!/bin/bash

#SBATCH --nodes 1
#SBATCH --job-name=hex_lattice
#SBATCH --time=2-12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
srun ./tolerance.sh --oversubscribe
