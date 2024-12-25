#!/bin/bash
#SBATCH --job-name=find_conda_path
#SBATCH --ntasks=1
#SBATCH --time=0-00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

# Load the Anaconda module
module load lang/Anaconda3/2024.02-1

# Print the location of conda
which conda
