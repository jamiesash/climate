#!/bin/bash
#SBATCH --job-name=stl_anomaly_monthly
#SBATCH --ntasks=1
#SBATCH --partition=shared
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --error=error-%A_%a.err
#SBATCH --output=output-%A_%a.out
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --mail-user=jamesash@hawaii.edu

# Load the Anaconda module
module load lang/Anaconda3/2024.02-1

# Initialize conda (this is needed to make conda available)
eval "$(conda shell.bash hook)"
conda activate blooms

# Manually source the conda initialization script

# Change to the script directory
cd /home/jamesash/climate/code/pycode

# Run the Python script
python3 chla_stl_mon_proc.py

