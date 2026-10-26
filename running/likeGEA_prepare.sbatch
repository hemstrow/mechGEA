#!/bin/bash
#SBATCH --job-name=likeGEA_prep
#SBATCH --output=likeGEA_prep.%j.out
#SBATCH --error=likeGEA_prep.%j.err
#SBATCH --mem=5G
#SBATCH --time=00:30:00
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --mail-type=END,FAIL
#SBATCH  --mail-user=ericacnr@colostate.edu

source ~/.bashrc
conda activate R

beagle_path="/scratch/alpine/ericacnr@colostate.edu/BCRF_Historical/lc_pipe/angsd/snp_calling_global"
docs_path="/scratch/alpine/ericacnr@colostate.edu/BCRF_Historical/lc_pipe/docs"

awk 'NR>1 && $5=="current" {
    print "Ind" NR-1
}' "$docs_path/hBCRF.locopipe.norep.noPP.noLP.three_sets.txt" \
> current_beagle_inds.txt


zcat "$beagle_path/combined.subsetted.beagle.gz" | head -1 |
awk '{for(i=4; i<=NF; i+=3) print $i}' |
sort > beagle_file_inds.txt

# Intended IDs
sort current_beagle_inds.txt > intended_inds.txt

# In intended list but missing from BEAGLE file
comm -23 intended_inds.txt beagle_file_inds.txt > missing_inds.txt

if [[ -s missing_inds.txt ]]; then
    echo "ERROR: Missing individuals:"
    cat missing_inds.txt
    exit 1
fi

zcat "$beagle_path/combined.subsetted.beagle.gz" |
awk '
BEGIN {
    while ((getline x < "current_beagle_inds.txt") > 0) {
        gsub(/\r/, "", x)
        keep[x] = 1
    }
}

NR == 1 {
    use[1] = use[2] = use[3] = 1

    for (i = 4; i <= NF; i += 3) {
        if ($i in keep) {
            use[i]   = 1
            use[i+1] = 1
            use[i+2] = 1
        }
    }
}

{
    first = 1
    for (i = 1; i <= NF; i++) {
        if (use[i]) {
            if (!first) printf "\t"
            printf "%s", $i
            first = 0
        }
    }
    printf "\n"
}
' | gzip > "$beagle_path/current.combined.subsetted.norep.unrel.noPP.noLP.beagle.gz"

#------------------
# Create chunks using likeGEA

echo "Running get_chunks()..."

Rscript likeGEA_chunks.R


#------------------
# Count chunks

chunk_dir="$beagle_path/likeGEA_chunks_3"

nchunks=$(find "$chunk_dir" -maxdepth 1 \
    -name 'chunk_*.beagle.gz' | wc -l)

if [[ "$nchunks" -eq 0 ]]; then
    echo "ERROR: No chunks created"
    exit 1
fi

echo "Created $nchunks chunks."

last_index=$((nchunks - 1))

#------------------
# Submit migration job

echo "Submitting migration job..."

sbatch \
    --export=ALL,LAST_INDEX="$last_index" \
    likeGEA_migration.sbatch
    
#------------------
#### OLD STUFF
# Submitting each chunk as it's own job
#chunk_dir="$beagle_path/likeGEA_chunks"

#nchunks=$(find "$chunk_dir" -maxdepth 1 \
    -name 'chunk_*.beagle.gz' | wc -l)

#if [[ "$nchunks" -eq 0 ]]; then
    echo "ERROR: No chunks created"
    exit 1
#fi

#echo "Created $nchunks chunks."

#last_index=$((nchunks - 1))

#echo "Submitting array 0-${last_index}%4"

#sbatch --array="0-${last_index}%4" likeGEA_array.sbatch

