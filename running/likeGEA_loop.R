library(likeGEA)
library(tictoc)
library(cmdstanr)
library(ggplot2)
library(data.table)
library(dplyr)

task_id <- as.integer(commandArgs(trailingOnly = TRUE)[1])
task_id <- 0

base <- "/scratch/alpine/ericacnr@colostate.edu/BCRF_Historical/lc_pipe/angsd/snp_calling_global"

chunk_dir <- file.path(base, "likeGEA_chunks")
result_dir <- file.path("/scratch/alpine/ericacnr@colostate.edu/BCRF_Historical/GEA/likeGEA_results")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

result_file <- file.path(
  result_dir,
  sprintf("likeGEA_chunk_%03d.rds", task_id)
)

if (file.exists(result_file)) {
  message("Already completed: ", result_file)
  quit(save = "no", status = 0)
}

chunk_file <- file.path(
  chunk_dir,
  sprintf("chunk_%03d.beagle.gz", task_id)
)

stopifnot(file.exists(chunk_file))

message(print(chunk_file))

likelihoods <- beagle2pl(chunk_file)

message(
  "Loaded ", nrow(likelihoods),
  " loci and ", ncol(likelihoods),
  " individuals"
)

indv <- fread("/scratch/alpine/ericacnr@colostate.edu/BCRF_Historical/lc_pipe/docs/hBCRF.locopipe.norep.noPP.noLP.three_sets.txt")

head(example_gl$patch_info)
patch_info <- indv %>% filter(set=="current") %>% select(sample_name, Population) %>% rename(ind_id=sample_name, patch=Population)

message("Beagle individuals:", head(colnames(likelihoods)))

message("Metadata individuals:", head(patch_info$ind_id))

meta <- read.csv("/scratch/alpine/ericacnr@colostate.edu/BCRF_Historical/GEA/BCRF_full_meta_table_propnames.csv") %>% rename(ind_id=sample_name)
patches_xy <- patch_info %>% left_join(meta %>% select(ind_id, Lat, Long), by="ind_id") %>% rename(lat=Lat, lon=Long) %>% select(-ind_id) %>%
  distinct(patch, .keep_all = TRUE)

env_pre <- read.csv("/scratch/alpine/ericacnr@colostate.edu/BCRF_Historical/GEA/current_climateNA_breeding_summary_scaled_matched_locations.csv")

env <- env_pre %>% select(-YEAR, -climate_period) %>% rename(patch=Population) %>% select(patch, PAS_tot)

stopifnot(
  !anyDuplicated(env$patch),
  setequal(env$patch, patches_xy$patch),
  setequal(patch_info$patch, patches_xy$patch)
)

# Arrange patches consistently
patch_info <- indv %>%
  mutate(beagle_id = paste0("Ind", row_number())) %>%
  filter(set == "current") %>%
  transmute(
    ind_id = beagle_id,
    patch = Population
  )

patches_xy <- patches_xy %>%
  arrange(patch)

env <- env %>%
  arrange(patch)

stopifnot(identical(env$patch, patches_xy$patch))

env_vars <- env %>%
  dplyr::select(-patch)

migration_loci <- 1000

task_dir <- file.path(
  result_dir,
  sprintf("work_chunk_%03d", task_id)
)

dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
setwd(task_dir)

res <- likeGEA(
  likelihoods,
  migration_loci = migration_loci,
  patch_info = patch_info,
  patches_xy = patches_xy,
  env = env
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

temp_file <- paste0(result_file, ".tmp")

saveRDS(res, temp_file)

if (!file.rename(temp_file, result_file)) {
  stop("Failed to save ", result_file)
}

message("Saved result: ", result_file)
message("Chunk ", task_id, " complete.")