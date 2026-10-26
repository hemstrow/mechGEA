# likeGEA using a PRECOMPUTED migration matrix, one genomic chunk at a time.
# Requires existing likeGEA helpers and Stan assets (see notes below).
# IMPORTANT: run each simultaneous SLURM array job in a separate working
# directory: the downstream Stan functions write ./stan and ./files files.

likeGEA_no_mig <- function(Phred, patch_info, patches_xy, env, migration,
                           migration_patches = NULL,
                           k_sd = 4, p0 = 20, seed = 42,
                           chains_gea = 2, refresh_gea = 200,
                           iter_warmup_gea = 250, iter_sampling_gea = 500,
                           plotting_CI = 0.95,
                           run_plot = TRUE, run_neff = TRUE,
                           par = max(1L, parallel::detectCores() - 2L),
                           verbose = TRUE) {
  # 1. Prepare this chunk's GEA-only data; no migration simulation or fit.
  stan_data <- prepare_gea_data(Phred, patch_info, patches_xy, env,
                                verbose = verbose)
  patch_names <- attr(stan_data, "patch_names")
  
  # 2. Check/align the existing migration matrix to the Stan patch order.
  Migration <- as.matrix(migration)
  if (!is.numeric(Migration) || nrow(Migration) != length(patch_names) ||
      ncol(Migration) != length(patch_names) || any(!is.finite(Migration)))
    stop("migration must be a finite numeric N_patch x N_patch matrix.")
  
  if (!is.null(rownames(Migration)) && !is.null(colnames(Migration))) {
    if (!setequal(rownames(Migration), patch_names) ||
        !setequal(colnames(Migration), patch_names))
      stop("Migration matrix patch names do not match patch_info.")
    Migration <- Migration[patch_names, patch_names, drop = FALSE]
  } else {
    # The original standalone migration code saves an unlabelled matrix.
    # You must supply the exact order used for the migration estimate.
    if (is.null(migration_patches) ||
        length(migration_patches) != length(patch_names) ||
        anyDuplicated(migration_patches) ||
        !setequal(migration_patches, patch_names))
      stop("Unlabelled migration matrix: provide migration_patches in its original row/column order.")
    dimnames(Migration) <- list(migration_patches, migration_patches)
    Migration <- Migration[patch_names, patch_names, drop = FALSE]
  }
  
  # 3. Run Stan GEA with the precomputed migration matrix.
  GEA_res <- run_stan_GEA_model(
    stan_data = stan_data, migration_results = Migration,
    k_sd = k_sd, p0 = p0, chains = chains_gea,
    refresh = refresh_gea, iter_warmup = iter_warmup_gea,
    iter_sampling = iter_sampling_gea, seed = seed, par = par
  )
  
  # 4. Plot and estimate N_eff if requested.
  if (run_plot) {
    plot_gea(stan_data = stan_data, migration_results = Migration,
             beta_draw = GEA_res$beta_draw,
             beta_summary = GEA_res$beta_summary,
             rho_kappa_mu = GEA_res$rho_kappa_mu,
             env_names = colnames(env)[-1], CI = plotting_CI)
  }
  neff <- if (run_neff) estimate_neff(stan_data) else NULL
  
  return(list(GEA_results = GEA_res,
              migration_results = Migration,
              Neff_results = neff,
              stan_data = stan_data))
}


prepare_gea_data <- function(Phred, patch_info, patches_xy, env,
                             verbose = TRUE) {
  # Matches the patch grouping convention used in migration_standalone.R.
  Phred <- as.matrix(Phred)
  required <- c("ind_id", "patch")
  if (!all(required %in% names(patch_info)))
    stop("patch_info must have ind_id and patch columns.")
  if (!all(c("patch", "lon", "lat") %in% names(patches_xy)))
    stop("patches_xy must have patch, lon, lat columns.")
  if (!"patch" %in% names(env)) stop("env must have a patch column first.")
  if (names(env)[1] != "patch")
    stop("The first column of env must be 'patch'.")
  if (ncol(Phred) != nrow(patch_info))
    stop("Phred columns must match patch_info rows.")
  if (anyNA(patch_info[, required]) || anyDuplicated(patch_info$ind_id))
    stop("Invalid individual IDs or patch assignments.")
  if (!is.null(colnames(Phred)) &&
      !identical(colnames(Phred), as.character(patch_info$ind_id)))
    stop("Phred column names are not in patch_info$ind_id order.")
  
  ord <- order(as.character(patch_info$patch))
  patch_info <- patch_info[ord, , drop = FALSE]
  Phred <- Phred[, ord, drop = FALSE]
  patches <- unique(as.character(patch_info$patch))
  N_table <- as.integer(table(factor(patch_info$patch, levels = patches)))
  N_patch <- length(patches)
  if (N_patch < 2L) stop("At least two patches are required.")
  
  # Ensure environments represent exactly the same patches, in the same order.