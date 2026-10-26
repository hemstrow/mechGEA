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
    stop("migration_loci exceeds loci available in this chunk.")
  if (N_sim < 1 || N_sim != as.integer(N_sim))
    stop("N_sim must be a positive integer.")
  
  set.seed(seed)
  migration_idx <- sort(sample.int(N_locus, size = migration_loci))
  Phred_mig <- Phred[migration_idx, , drop = FALSE]
  rm(Phred)
  
  # Only convert the chosen loci, never the entire input chunk.
  values <- as.vector(Phred_mig)
  valid <- !is.na(values) & stringi::stri_detect_regex(values, "^[0-9]+,[0-9]+,[0-9]+$")
  values[!valid] <- "."
  Phred_mig[] <- values
  rm(values, valid)
  if (verbose) {
    cat(sprintf("Selected %d of %d loci across %d patches.\n",
                migration_loci, N_locus, N_patch))
    cat("Converting selected Phred scores to genotype likelihoods ...\n")
  }
  gl_migration <- likeGEA:::phred_to_loglik(Phred_mig)
  rm(Phred_mig)
  
  set.seed(seed + 1L)
  pi_off_sum <- matrix(0, N_patch, N_patch)
  pi_diag_sum <- numeric(N_patch)
  for (s in seq_len(N_sim)) {
    geno_s <- likeGEA:::genotype_from_gl(gl_migration)
    pi_s <- likeGEA:::compute_pi_genotype(geno_s, N_table)
    pi_off_sum <- pi_off_sum + pi_s$pi_off
    pi_diag_sum <- pi_diag_sum + pi_s$pi_diag
    if (verbose && (s %% 10L == 0L || s == N_sim))
      cat(sprintf("  Resample %d / %d\n", s, N_sim))
  }
  
  # geodist expects longitude first, latitude second; returns meters.
  DistM <- geodist::geodist(coords[, c("lon", "lat")], measure = "geodesic")
  if (any(!is.finite(DistM))) stop("Distance matrix has non-finite values.")
  if (sd(DistM[upper.tri(DistM)]) == 0)
    stop("All pairwise geographic distances are identical; cannot standardize.")
  
  migration_data <- list(
    N_patch = N_patch,
    patches = patches,
    N_locus_neutral = migration_loci,
    N_table = N_table,
    DistM = DistM,
    pi_off = pi_off_sum / N_sim,
    pi_diag = pi_diag_sum / N_sim
  )
  dir.create("files", showWarnings = FALSE, recursive = TRUE)
  saveRDS(migration_data, "files/1_migration_data.RDS")
  if (verbose) cat("Saved files/1_migration_data.RDS\n")
  return(migration_data)
}


estimate_migration <- function(data,
                               model = "migration_fast",
                               max_attempts = 10,
                               target_div_pct = 1,
                               target_rhat_pct = 1,
                               chains = 4,
                               iter_sampling = 500,
                               iter_warmup = 250,
                               max_treedepth = 12,
                               refresh = 200,
                               cores = max(1L, parallel::detectCores() - 2L)) {
  
  # Copy the original likeGEA Stan files into ./stan before running.
  stan_dir <- "/projects/ericacnr@colostate.edu/likeGEA_bebug/likeGEA/stan/"
  
  stan_file <- file.path(stan_dir, paste0(sub("\\.stan$", "", model), ".stan"))
  hpp_file <- file.path(stan_dir, "coalescent_time.hpp")
  
  stopifnot(file.exists(stan_file), file.exists(hpp_file))
  
  N_patch <- data$N_patch
  stopifnot(is.matrix(data$pi_off),
            all(dim(data$pi_off) == c(N_patch, N_patch)),
            length(data$pi_diag) == N_patch,
            all(is.finite(data$pi_off[upper.tri(data$pi_off)])),
            all(is.finite(data$pi_diag)))
  
  scale_sd <- sd(data$DistM[upper.tri(data$DistM)])
  if (!is.finite(scale_sd) || scale_sd <= 0)
    stop("Cannot standardize geographic distances.")
  DistM_std <- data$DistM / scale_sd
  covariates <- matrix(likeGEA:::convert_T_to_t(DistM_std), ncol = 1)
  if (any(!is.finite(covariates)))
    stop("Migration covariates contain missing or non-finite values.")
  
  mod <- cmdstanr::cmdstan_model(
    stan_file,
    cpp_options = list(
      USER_HEADER = normalizePath(hpp_file),
      stan_threads = TRUE
    ),
    stanc_options = list(
      "allow-undefined" = TRUE,
      "O1"
    )
  )
  
  threads_per_chain <- max(1L, floor(cores / chains))
  dat <- list(
    pi_off = likeGEA:::convert_T_to_t(data$pi_off),
    pi_diag = data$pi_diag,
    covariates = covariates,
    N_cov = 1L,
    N_patch = N_patch,
    gamma_sd = 1
  )
  
  total_iter <- chains * iter_sampling
  best_fit <- NULL
  best_score <- Inf
  best_div_pct <- NA_real_
  best_rhat_pct <- NA_real_
  best_seed <- NA_integer_
  
  for (attempt in seq_len(max_attempts)) {
    fit_seed <- as.integer(Sys.time()) + attempt
    cat(sprintf("\n=== %s attempt %d/%d (seed=%d) ===\n",
                model, attempt, max_attempts, fit_seed))
    fit <- mod$sample(
      data = dat, chains = chains, parallel_chains = chains,
      seed = fit_seed,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      threads_per_chain = threads_per_chain,
      max_treedepth = max_treedepth,
      refresh = refresh,
      show_messages = TRUE,
      show_exceptions = FALSE
    )
    
    n_div <- sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
    div_pct <- 100 * n_div / total_iter
    rhat_vec <- fit$summary()$rhat
    rhat_vec <- rhat_vec[is.finite(rhat_vec)]
    rhat_pct <- if (length(rhat_vec)) 100 * mean(rhat_vec > 1.05) else Inf
    cat(sprintf("Divergences: %.2f%%; Rhat > 1.05: %.2f%%\n",
                div_pct, rhat_pct))
    
    score <- max(div_pct, rhat_pct)
    if (score < best_score) {
      best_fit <- fit
      best_score <- score
      best_div_pct <- div_pct
      best_rhat_pct <- rhat_pct
      best_seed <- fit_seed
    }
    if (div_pct < target_div_pct && rhat_pct < target_rhat_pct) break
  }
  if (is.null(best_fit)) stop("No successful migration fit was obtained.")
  if (best_div_pct >= target_div_pct || best_rhat_pct >= target_rhat_pct)
    warning(sprintf("Targets not reached; using best fit (seed=%d, divergences=%.2f%%, high Rhat=%.2f%%).",
                    best_seed, best_div_pct, best_rhat_pct))
  
  M_summary <- best_fit$summary(variables = "Migration")
  Migration <- matrix(M_summary$mean, nrow = N_patch, ncol = N_patch)
  if (!is.null(data$patches))
    dimnames(Migration) <- list(data$patches, data$patches)
  saveRDS(Migration, "files/2_migration.RDS")
  cat("Saved files/2_migration.RDS\n")
  cat(sprintf("Migration: max Rhat=%.3f, min bulk ESS=%.0f\n",
              max(M_summary$rhat, na.rm = TRUE),
              min(M_summary$ess_bulk, na.rm = TRUE)))
  return(Migration)
}
