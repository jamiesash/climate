#!/bin/bash
#SBATCH --job-name=anom_chl
#SBATCH --ntasks=1
#SBATCH --partition=shared
#SBATCH --time=1-12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=300G
#SBATCH --error=error-%A_%a.err
#SBATCH --output=output-%A_%a.out
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --mail-user=jamesash@hawaii.edu

# Load the Anaconda module
module load lang/Anaconda3/2024.02-1

# Initialize conda (this is needed to make conda available)
eval "$(conda shell.bash hook)"

# This env contains a user space installation of R.
# Make sure to not use the R module as that will clash. 
conda activate rblooms

# Change to the script directory
cd /home/jamesash/climate/code/rcode/anomaly

# Run the Python script
Rscript chla_proc_day.R

