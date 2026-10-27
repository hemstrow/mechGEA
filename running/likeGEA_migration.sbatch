#!/bin/bash
#SBATCH --job-name=likeGEA_mig
#SBATCH --output=likeGEA_mig.%j.out
#SBATCH --error=likeGEA_mig.%j.err
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ericacnr@colostate.edu

source ~/.bashrc
conda activate R

migration_file="/scratch/alpine/ericacnr@colostate.edu/BCRF_Historical/GEA/likeGEA_results/migration/Migration_shared.RDS"
rm -f "$migration_file"

echo "Running migration estimation..."

Rscript likeGEA_migration.R

#------------------
# Check migration completed successfully

if [[ ! -s "$migration_file" ]]; then
    echo "ERROR: Migration_shared.RDS does not exist or is empty."
    echo "Migration failed. GEA array will NOT be submitted."
    exit 1
fi

echo "Migration output found: $migration_file"

# #------------------
# # Submit GEA array
# 
# echo "Submitting GEA array 0-${LAST_INDEX}%4"
# 
# sbatch \
#     --array="0-${LAST_INDEX}%4" \
#     likeGEA_array.sbatch