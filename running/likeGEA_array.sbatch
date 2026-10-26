#!/bin/bash
#SBATCH --job-name=likeGEA_run
#SBATCH --output=likeGEA_run.%A_%a.out
#SBATCH --error=likeGEA_run.%A_%a.err
#SBATCH --mem=46G
#SBATCH --time=23:59:00
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --mail-type=END,FAIL
#SBATCH  --mail-user=ericacnr@colostate.edu

source ~/.bashrc
conda activate R

#sbatch --array=0-7 likeGEA_array.sbatch

Rscript likeGEA_running.R "$SLURM_ARRAY_TASK_ID"