#!/usr/bin/env Rscript
# Run one round-robin Beagle chunk through GEA using the shared migration matrix.
# Usage: Rscript likeGEA_running.R 0   (or pass $SLURM_ARRAY_TASK_ID)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(likeGEA)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !grepl("^[0-9]+$", args[1]))
  stop("Supply one zero-based chunk index: Rscript likeGEA_running.R 0")
task_id <- as.integer(args[1])

# Source alongside this R script, regardless of SLURM's launch directory.
source("../R/likeGEA_no_mig.R")

base <- "/scratch/alpine/ericacnr@colostate.edu/BCRF_Historical/lc_pipe/angsd/snp_calling_global"
root <- "/scratch/alpine/ericacnr@colostate.edu/BCRF_Historical"
chunk_dir <- file.path(base, "likeGEA_chunks_3")
#representative_draw/ for the other version
migration_dir <- file.path(root, "GEA/likeGEA_results_3/migration/")
result_dir <- file.path(root, "GEA/likeGEA_results_sensativity/posterior_representative")
migration_file <- file.path(migration_dir, "Migration_representative.RD")
chunk_file <- file.path(chunk_dir, sprintf("chunk_%03d.beagle.gz", task_id))
result_file <- file.path(result_dir, sprintf("likeGEA_chunk_%03d.rds", task_id))
task_dir <- file.path(result_dir, sprintf("work_chunk_%03d", task_id))

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(result_file)) {
  message("Already completed: ", result_file)
  quit(save = "no", status = 0)
}
if (!file.exists(chunk_file)) stop("Chunk not found: ", chunk_file)
if (!file.exists(migration_file)) stop("Shared migration matrix not found: ", migration_file)

# Import the exact individual IDs encoded in the filtered Beagle file.
message("Loading ", chunk_file)

header <- strsplit(
  readLines(gzfile(chunk_file), n = 1),
  "\\s+"
)[[1]]

# Beagle has 3 columns per individual
beagle_ids <- header[seq(4, length(header), by = 3)]

# Convert genotype likelihoods
likelihoods <- beagle2pl(chunk_file)

# Restore individual IDs
stopifnot(ncol(likelihoods) == length(beagle_ids))

colnames(likelihoods) <- beagle_ids

message("Loaded ", nrow(likelihoods), " loci and ",
        ncol(likelihoods), " individuals")

indv <- as.data.frame(
  fread(file.path(root, "lc_pipe/docs/hBCRF.locopipe.norep.noPP.noLP.three_sets.txt"))
)

# IMPORTANT: Ind IDs count every row in the ORIGINAL metadata, before filtering.
patch_info <- indv %>%
  mutate(beagle_id = paste0("Ind", row_number())) %>%
  filter(set == "current") %>%
  transmute(ind_id = beagle_id, patch = as.character(Population)) %>%
  arrange(patch)

stopifnot(!anyDuplicated(patch_info$ind_id), !anyNA(patch_info$patch))
if (is.null(colnames(likelihoods)))
  stop("Beagle-derived Phred matrix must retain its individual column names.")
if (!setequal(colnames(likelihoods), patch_info$ind_id)) {
  stop("Beagle IDs do not match current metadata IDs. Check header, filtering and Ind numbering.")
}
likelihoods <- likelihoods[, patch_info$ind_id, drop = FALSE]
stopifnot(identical(colnames(likelihoods), patch_info$ind_id))

# Geographic coordinates from the original sample IDs, averaged within patch.
meta <- read.csv(file.path(root, "GEA/BCRF_full_meta_table_propnames.csv")) %>%
  rename(sample_name = BGP_ID)
xy_ind <- indv %>%
  filter(set == "current") %>%
  select(sample_name, Population) %>%
  rename(patch = Population) %>%
  left_join(
    meta %>% select(sample_name, Lat, Long),
    by = "sample_name"
  )
if (anyNA(xy_ind$Lat) || anyNA(xy_ind$Long))
  stop("Missing sample coordinates; check metadata join before using patch means.")
patches_xy <- xy_ind %>%
  group_by(patch) %>%
  summarise(lon = mean(Long), lat = mean(Lat), .groups = "drop") %>%
  arrange(patch)

# Keep the same environmental predictor used in the original chunk pipeline.
env_pre <- read.csv(file.path(root, "GEA/current_climateNA_breeding_summary_scaled_matched_locations.csv"))
env <- env_pre %>%
  select(-YEAR, -climate_period) %>%
  rename(patch = Population) %>%
  select(patch, PAS_tot) %>%
  arrange(patch)

stopifnot(!anyDuplicated(env$patch),
          identical(unique(patch_info$patch), patches_xy$patch),
          identical(env$patch, patches_xy$patch),
          all(is.finite(patches_xy$lon)),
          all(is.finite(patches_xy$lat)),
          all(is.finite(env$PAS_tot)))

Migration <- readRDS(migration_file)
if (is.null(rownames(Migration)) || is.null(colnames(Migration)))
  stop("The shared migration matrix must have patch row/column names. Re-save the named matrix.")
if (!setequal(rownames(Migration), patches_xy$patch) ||
    !setequal(colnames(Migration), patches_xy$patch))
  stop("Migration matrix patches do not match this chunk's patches.")

# Separate work directory is essential: all jobs write files/ and stan/.
dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
setwd(task_dir)
message("Working in ", getwd())

cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
if (is.na(cores) || cores < 1L) cores <- 1L

res <- likeGEA_no_mig(
  Phred = likelihoods,
  patch_info = patch_info,
  patches_xy = patches_xy,
  env = env,
  migration = Migration,
  chains_gea = 2,
  par = cores,
  run_plot = TRUE,
  run_neff = TRUE
)

# Write a completed result only after all requested stages succeed.
tmp_file <- paste0(result_file, ".", Sys.getpid(), ".tmp")
saveRDS(res, tmp_file)
if (!file.rename(tmp_file, result_file))
  stop("Could not finalize result file: ", result_file)
message("Completed and saved: ", result_file)
