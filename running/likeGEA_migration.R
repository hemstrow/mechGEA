
library(data.table)
library(dplyr)
library(likeGEA)

source("migration_standalone_updated.R")

base <- "/scratch/alpine/ericacnr@colostate.edu/BCRF_Historical/lc_pipe/angsd/snp_calling_global"

result_dir <- "/scratch/alpine/ericacnr@colostate.edu/BCRF_Historical/GEA/likeGEA_results_3"

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

#------------------
# Load first round-robin chunk

chunk_file <- file.path(
  base, "likeGEA_chunks_3", "chunk_000.beagle.gz"
)

stopifnot(file.exists(chunk_file))

likelihoods <- beagle2pl(chunk_file)

message("Loaded ", nrow(likelihoods), " loci and ",
        ncol(likelihoods), " individuals")

#------------------
# Prepare patch_info

indv <- fread(
  "/scratch/alpine/ericacnr@colostate.edu/BCRF_Historical/lc_pipe/docs/hBCRF.locopipe.norep.noPP.noLP.three_sets.txt"
)

patch_info <- indv %>%
  mutate(beagle_id = paste0("Ind", row_number())) %>%
  filter(set == "current") %>%
  transmute(
    ind_id = beagle_id,
    patch = Population
  )

# Match metadata to Beagle columns
stopifnot(
  !anyDuplicated(patch_info$ind_id),
  setequal(colnames(likelihoods), patch_info$ind_id)
)

# Arrange individuals by patch
patch_info <- patch_info %>%
  arrange(patch)

likelihoods <- likelihoods[
  , patch_info$ind_id, drop = FALSE
]

stopifnot(identical(colnames(likelihoods), patch_info$ind_id))

#------------------
# Prepare patches_xy

meta <- read.csv(
  "/scratch/alpine/ericacnr@colostate.edu/BCRF_Historical/GEA/BCRF_full_meta_table_propnames.csv"
) %>%
  rename(ind_id = sample_name)

patches_xy <- indv %>%
  filter(set == "current") %>%
  select(sample_name, Population) %>%
  rename(ind_id = sample_name, patch = Population) %>%
  left_join(
    meta %>% select(ind_id, Lat, Long),
    by = "ind_id"
  ) %>%
  group_by(patch) %>%
  summarise(
    lon = mean(Long, na.rm = TRUE),
    lat = mean(Lat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(patch)

stopifnot(
  identical(unique(patch_info$patch), patches_xy$patch),
  all(is.finite(patches_xy$lon)),
  all(is.finite(patches_xy$lat))
)

#------------------
# Run migration

migration_dir <- file.path(result_dir, "migration")
dir.create(migration_dir, recursive = TRUE, showWarnings = FALSE)

setwd(migration_dir)

Migration <- migration(
  Phred = likelihoods,
  patch_info = patch_info,
  patches_xy = patches_xy,
  migration_loci = 50000,
  N_sim_mig = 100,
  chains_mig = 4,
  cores_mig = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
)

# Explicitly name the matrix for future GEA chunks
dimnames(Migration) <- list(
  patches_xy$patch,
  patches_xy$patch
)

saveRDS(Migration, "Migration_shared.RDS")

message("Migration estimation complete.")
