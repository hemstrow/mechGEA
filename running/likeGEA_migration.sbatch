#!/bin/bash
#SBATCH --job-name=likeGEA_mig
#SBATCH --output=likeGEA_mig.%j.out
#SBATCH --error=likeGEA_mig.%j.err
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ericacnr@colostate.edu

set -euo pipefail

source ~/.bashrc
conda activate R

echo "Running migration estimation..."

Rscript likeGEA_migration.R

echo "Migration completed successfully."

#------------------
# Submit GEA array

echo "Submitting GEA array 0-${LAST_INDEX}%4"

sbatch \
    --array="0-${LAST_INDEX}%4" \
    likeGEA_array.sbatch