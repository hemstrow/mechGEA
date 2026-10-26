# Standalone migration workflow from a genome-wide round-robin Beagle chunk.
# Required existing likeGEA helper functions (source or load them):
#   phred_to_loglik(), genotype_from_gl(), compute_pi_genotype(),
#   convert_T_to_t().
# Required Stan assets: stan/migration_fast.stan (or migration_M3b.stan)
# and stan/coalescent_time.hpp, relative to the R working directory.
# Required packages: stringi, geodist, cmdstanr.

migration <- function(Phred, patch_info, patches_xy,
                      migration_loci = 50000,
                      N_sim_mig = 100,
                      seed = 42,
                      verbose = TRUE,
                      migration_model = "migration_fast",
                      max_attempts_mig = 10,
                      target_div_pct_mig = 1,
                      target_rhat_pct_mig = 1,
                      chains_mig = 4,
                      iter_sampling_mig = 500,
                      iter_warmup_mig = 250,
                      max_treedepth_mig = 12,
                      refresh_mig = 200,
                      cores_mig = max(1L, parallel::detectCores() - 2L)) {
  
  mig_data <- prepare_mig_data(
    Phred = Phred,
    individuals = patch_info,
    coordinates = patches_xy,
    migration_loci = migration_loci,
    N_sim = N_sim_mig,
    seed = seed,
    verbose = verbose
  )
  
  Migration <- estimate_migration(
    data = mig_data,
    model = migration_model,
    max_attempts = max_attempts_mig,
    target_div_pct = target_div_pct_mig,
    target_rhat_pct = target_rhat_pct_mig,
    chains = chains_mig,
    iter_sampling = iter_sampling_mig,
    iter_warmup = iter_warmup_mig,
    max_treedepth = max_treedepth_mig,
    refresh = refresh_mig,
    cores = cores_mig
  )
  
  return(Migration)
}


prepare_mig_data <- function(Phred, individuals, coordinates,
                             migration_loci = 50000, N_sim = 100,
                             seed = 42, verbose = TRUE) {
  
  # Phred: loci x individuals, with cells like "12,34,56" or ".".
  # individuals: data frame with ind_id and patch, in Phred column order.
  # coordinates: one row per patch, with patch, lon, lat.
  Phred <- as.matrix(Phred)
  stopifnot(all(c("ind_id", "patch") %in% names(individuals)),
            all(c("patch", "lon", "lat") %in% names(coordinates)))
  if (ncol(Phred) != nrow(individuals))
    stop("Phred column count must equal number of rows in individuals.")
  if (anyNA(individuals[, c("ind_id", "patch")]) ||
      anyDuplicated(individuals$ind_id))
    stop("Individual IDs and patch assignments must be nonmissing; IDs unique.")
  if (!is.null(colnames(Phred)) &&
      !identical(colnames(Phred), as.character(individuals$ind_id)))
    stop("Phred column names do not match individuals$ind_id order.")
  
  # Group individuals by patch: required by compute_pi_genotype().
  # Reorder Phred columns together with metadata to maintain alignment.
  ord <- order(as.character(individuals$patch))
  individuals <- individuals[ord, , drop = FALSE]
  Phred <- Phred[, ord, drop = FALSE]
  patches <- unique(as.character(individuals$patch))
  N_table <- as.integer(table(factor(individuals$patch, levels = patches)))
  N_patch <- length(patches)
  if (N_patch < 2L) stop("At least two patches are required.")
  
  coords <- as.data.frame(coordinates)
  if (anyDuplicated(coords$patch) ||
      !setequal(as.character(coords$patch), patches))
    stop("coordinates$patch must contain every sampled patch exactly once.")
  coords <- coords[match(patches, as.character(coords$patch)), , drop = FALSE]
  if (any(!is.finite(as.numeric(coords$lon))) ||
      any(!is.finite(as.numeric(coords$lat))))
    stop("Coordinates contain missing or non-finite lon/lat values.")
  
  N_locus <- nrow(Phred)
  if (N_locus == 0L) stop("Phred has no loci.")
  if (length(migration_loci) != 1L || !is.finite(migration_loci) ||
      migration_loci <= 0 || migration_loci != as.integer(migration_loci))
    stop("migration_loci must be a positive integer.")
  if (migration_loci > N_locus)