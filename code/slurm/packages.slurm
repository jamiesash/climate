#!/bin/bash
#SBATCH --job-name=cli_map

#SBATCH --ntasks=1

#SBATCH --partition=shared
## 3 day max run time for public partitions, except 4 hour limit  in sandbox
#SBATCH --time=2-00:00:00 ## time format is DD-HH:MM:SS

#SBATCH --cpus-per-task=8
#SBATCH --mem=1000G ## max amount of memory per node you require

#SBATCH --error=hello-%A_%a.err ## %A - filled with jobid
#SBATCH --output=hello-%A_%a.out ## %A - filled with jobid

##SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
##SBATCH --mail-user=user@test.org
##SBATCH --array=1-10
##SBATCH --core-spec=0 ## Uncomment to allow jobs to request all cores on a node    

## All options and environment variables found on schedMD site: http://slurm.schedmd.com/sbatch.html

## module purge
module load lang/Anaconda3/2023.03-1
source activate blooms

Rscript libraries_setup.R
