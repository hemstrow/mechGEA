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
  if (anyNA(env$patch) || anyDuplicated(env$patch) ||
      !setequal(as.character(env$patch), patches))
    stop("env must contain exactly one row for each sampled patch.")
  if (anyNA(patches_xy$patch) || anyDuplicated(patches_xy$patch) ||
      !setequal(as.character(patches_xy$patch), patches))
    stop("patches_xy must contain exactly one row per sampled patch.")
  env <- env[match(patches, as.character(env$patch)), , drop = FALSE]
  Env <- as.matrix(env[, -1, drop = FALSE])
  if (!is.numeric(Env) || !ncol(Env) || any(!is.finite(Env)))
    stop("Environmental covariates must be numeric, complete, and finite.")
  N_locus <- nrow(Phred)
  if (!N_locus) stop("Phred contains zero loci.")
  
  if (verbose) cat(sprintf("Preparing %d loci from %d individuals in %d patches.\n",
                           N_locus, ncol(Phred), N_patch))
  values <- as.vector(Phred)
  valid <- !is.na(values) & stringi::stri_detect_regex(
    values, "^[0-9]+,[0-9]+,[0-9]+$")
  if (verbose) cat(sprintf("Missing or malformed Phred cells: %d / %d\n",
                           sum(!valid), length(valid)))
  values[!valid] <- "."
  Phred[] <- values
  rm(values, valid)
  
  gl_all <- likeGEA:::phred_to_loglik(Phred)
  rm(Phred)
  eff <- likeGEA:::compute_N_eff_G_eff_gl(gl_all, N_table = N_table)
  bb <- likeGEA:::compute_BetaBin_theta_gl(gl_all, N_table)
  G_eff <- eff$G_eff
  G_eff[is.na(G_eff)] <- 0
  
  stan_data <- list(
    N_patch = N_patch, N_locus = N_locus, N_env = ncol(Env),
    Env = Env, G_call = G_eff, TwoN = eff$TwoN_eff,
    BetaBin_theta = bb$BetaBin_theta, Group = bb$Group,
    N_table = N_table
  )
  attr(stan_data, "patch_names") <- patches
  dir.create("files", recursive = TRUE, showWarnings = FALSE)
  saveRDS(stan_data, "files/1_stan_data.RDS")
  if (verbose) cat("Saved files/1_stan_data.RDS\n")
  return(stan_data)
}


#' Run the core STAN model to estimate genotype-environment association with
#' likelihoods.
#'
#' Consider running via \code{\link{likeGEA}} for the full pipeline. Generally
#' for internal use or when partial re-runs are needed.
#'
#' @export
run_stan_GEA_model <- function(stan_data,
                               migration_results,
                               k_sd = 4,
                               p0 = 20,
                               refresh = 100,
                               chains = 2,
                               iter_warmup = 250,
                               iter_sampling = 500,
                               seed = 42,
                               par = parallel::detectCores() - 2){
  tictoc::tic("estimating Stan time")
  
  ##   migration_fast.stan -- Petkova-centered combination (default)
  ##   migration_M3b.stan  -- joint (pi_off, pi_diag) fit with diagonal anchor
  .source_model <- getFromNamespace(".source_model", "likeGEA")
  STAN_MODEL <- .source_model("structured.stan")
  if(!dir.exists("./stan")){dir.create("stan")}
  
  if(!file.exists(file.path("stan", "structured.stan"))){
    file.copy(STAN_MODEL, "./stan/structured.stan")
  }
  
  STAN_MODEL <- "stan/structured.stan"
  ## ---- Load data ------------------------------------------------------------
  
  data      <- stan_data; rm(stan_data)
  Migration <- migration_results; rm(migration_results)
  
  ## ---- Append hyperparameters -----------------------------------------------
  
  data$k_sd <- k_sd
  data$p0   <- p0
  
  ## structured.stan declares G_call and TwoN as real (lbeta + 1F1 arithmetic).
  data$G_call <- as.double(data$G_call)
  data$TwoN   <- as.double(data$TwoN)
  dim(data$G_call) <- c(data$N_locus, data$N_patch)
  dim(data$TwoN)   <- c(data$N_locus, data$N_patch)
  
  ## ---- Write data include file ----------------------------------------------
  
  likeGEA:::write_common_data_stan(
    as.matrix(data$Env),
    as.matrix(Migration),
    data$k_sd,
    "stan/structured_data.stan"
  )
  
  ## ---- Compile Stan model ---------------------------------------------------
  hpp <- .source_model("1F1_fast.hpp")
  nhpp <- file.path("stan", "1F1_fast.hpp")
  if(!file.exists(nhpp)){
    file.copy(hpp, nhpp)
  } else {
    cat("Using existing stan/1F1_fast.hpp\n")
  }
  hpp <- nhpp
  hpp_file <- normalizePath(hpp)
  
  mod <- cmdstanr::cmdstan_model(
    STAN_MODEL,
    cpp_options     = list(USER_HEADER = hpp_file, stan_threads = TRUE),
    stanc_options   = list("allow-undefined" = TRUE),
    include_paths   = normalizePath("stan"),
    force_recompile = FALSE
  )
  
  ## ---- Sample ---------------------------------------------------------------
  n_cores           <- par
  threads_per_chain <- floor(n_cores/chains)
  threads_per_chain <- max(c(threads_per_chain, 1))
  
  
  # DEV NOTE: should modify the function to take these as args at some point
  init_fun <- function() {
    list(
      mu          = runif(data$N_locus, 0.4, 0.6),
      kappa       = exp(rnorm(data$N_locus, log(1), 0.3)),
      rho_tilda   = exp(rnorm(data$N_locus, log(1), 0.2)),
      tau_unif    = runif(1, base::pi / 10, base::pi / 4),
      c_square    = abs(rnorm(1, 1.3, 0.3)),
      z_beta      = matrix(rnorm(data$N_env * data$N_locus, 0, 0.1),
                           nrow = data$N_env, ncol = data$N_locus),
      lambda_unif = matrix(runif(data$N_env * data$N_locus,
                                 base::pi / 10, base::pi / 5),
                           nrow = data$N_env, ncol = data$N_locus)
    )
  }
  
  fit <- mod$sample(
    data              = data,
    chains            = chains,
    seed              = seed,
    parallel_chains   = chains,
    init              = init_fun,
    show_exceptions   = TRUE,
    show_messages     = TRUE,
    refresh           = refresh,
    iter_warmup       = iter_warmup,
    threads_per_chain = threads_per_chain,
    iter_sampling     = iter_sampling
  )
  
  ## ---- Save posterior draws -------------------------------------------------
  
  dir.create("files", showWarnings = FALSE, recursive = TRUE)
  
  rho_kappa_mu <- fit$draws(variables = c("rho", "kappa", "mu"),
                            format = "draws_matrix")
  saveRDS(rho_kappa_mu, "files/4_rho_kappa_mu_draw.RDS")
  
  beta_draw <- fit$draws(variables = c("Beta"), format = "draws_matrix")
  saveRDS(beta_draw, "files/4_beta_draw.RDS")
  
  tau_draw <- fit$draws(variables = c("Tau"), format = "draws_matrix")
  saveRDS(tau_draw, "files/4_tau_draw.RDS")
  
  beta_summary <- likeGEA:::summarise_draws_parallel(beta_draw)
  saveRDS(beta_summary, "files/4_beta_summary.RDS")
  
  mlog_p <- likeGEA:::score_stan(beta_draw, data$N_locus, data$N_env)
  saveRDS(mlog_p, "files/4_mlog_p.RDS")
  
  cat("Saved files/4_{rho_kappa_mu_draw, beta_draw, tau_draw,",
      "beta_summary, mlog_p}.RDS\n")
  
  tictoc::toc(quiet = FALSE, func.toc = likeGEA:::toc_hours)
  
  return(list(rho_kappa_mu = rho_kappa_mu,
              beta_draw = beta_draw,
              tau = tau_draw,
              beta_summary = beta_summary,
              mlog_p = mlog_p))
  
}

#' Estimate effective sample sizes.
#'
#' Consider running via \code{\link{likeGEA}} for the full pipeline. Generally
#' for internal use or when partial re-runs are needed.
#'
#' @export
estimate_neff <- function(stan_data){
  tictoc::tic("estimating Neff time")
  
  ## ---- Load -----------------------------------------------------------------
  N_patch      <- stan_data$N_patch
  N_locus_tot  <- stan_data$N_locus
  G_eff_ana    <- stan_data$G_call         # analytical
  TwoN_eff_ana <- stan_data$TwoN           # analytical
  
  ## Subsample loci for speed (gl_p_est sampling is the cost driver).
  N_loci_use <- min(1000L, N_locus_tot)
  if (N_loci_use < N_locus_tot) {
    set.seed(789)
    loci_used <- sort(sample(N_locus_tot, N_loci_use))
    cat(sprintf("Subsampled %d of %d loci for N_eff estimation.\n",
                N_loci_use, N_locus_tot))
  } else {
    loci_used <- seq_len(N_locus_tot)
  }
  
  dat <- list(
    N_patch       = N_patch,
    N_locus       = length(loci_used),
    N_table       = as.integer(stan_data$N_table),
    BetaBin_theta = stan_data$BetaBin_theta[loci_used, , drop = FALSE],
    Group         = stan_data$Group        [loci_used, , drop = FALSE]
  )
  
  ## ---- Compile + sample gl_p_est.stan ---------------------------------------
  
  STAN_MODEL <- .source_model("gl_p_est.stan")
  if(!dir.exists("./stan")){dir.create("stan")}
  if(!file.exists("stan/gl_p_est.stan")){
    file.copy(STAN_MODEL, "./stan/gl_p_est.stan")
  }
  
  mod <- cmdstanr::cmdstan_model(
    "stan/gl_p_est.stan",
    cpp_options   = list(stan_threads = TRUE),
    stanc_options = list("O1")
  )
  
  n_cores           <- parallel::detectCores() - 2
  chain             <- 2
  threads_per_chain <- max(1L, floor(n_cores / chain))
  
  # DEV NOTE: Is it worth it to allow all of these to change?
  tictoc::tic("NUTS posterior on p")
  fit <- mod$sample(
    data              = dat,
    chains            = chain,
    parallel_chains   = chain,
    seed              = 546,
    threads_per_chain = threads_per_chain,
    refresh           = 100,
    iter_warmup       = 200,
    iter_sampling     = 500,
    show_messages     = TRUE
  )
  tictoc::toc(quiet = FALSE, func.toc = toc_hours)
  
  ## ---- Extract per-(locus, patch) posterior moments -------------------------
  
  p_summary <- summarise_draws_parallel(fit$draws("p", format = "draws_matrix"))
  idx_mat   <- do.call(rbind, regmatches(
    p_summary$variable,
    regexec("p\\[([0-9]+),([0-9]+)\\]", p_summary$variable)))
  i_idx     <- as.integer(idx_mat[, 2])
  k_idx     <- as.integer(idx_mat[, 3])
  stopifnot(length(i_idx) == length(loci_used) * N_patch,
            all(!is.na(i_idx)), all(!is.na(k_idx)))
  
  mu_post  <- matrix(NA_real_, length(loci_used), N_patch)
  sd_post  <- matrix(NA_real_, length(loci_used), N_patch)
  mu_post[cbind(i_idx, k_idx)] <- p_summary$mean
  sd_post[cbind(i_idx, k_idx)] <- p_summary$sd
  var_post <- sd_post^2
  
  ## ---- Beta(G+1, 2N-G+1) moment-match  (same convention as analytical) ------
  
  TwoN_eff_stan <- pmax(mu_post * (1 - mu_post) / pmax(var_post, 1e-20) - 3, 1e-6)
  G_eff_stan    <- (2 + TwoN_eff_stan) * mu_post - 1
  G_eff_stan    <- pmax(pmin(G_eff_stan, TwoN_eff_stan), 0)
  
  ## ---- Save -----------------------------------------------------------------
  
  
  
  ## Mean per-individual N_eff diagnostic (Stan estimate, like the original).
  Neff_mean <- mean(rowSums(TwoN_eff_stan) / 2) / sum(stan_data$N_table)
  cat(sprintf("Mean N_eff / N (Stan, across subsampled loci): %.4f\n", Neff_mean))
  
  ret <- list(
    loci_used     = loci_used,
    G_eff_stan    = G_eff_stan,
    TwoN_eff_stan = TwoN_eff_stan,
    G_eff_ana     = G_eff_ana   [loci_used, , drop = FALSE],
    TwoN_eff_ana  = TwoN_eff_ana[loci_used, , drop = FALSE],
    Neff_mean     = Neff_mean
  )
  
  dir.create("files", showWarnings = FALSE, recursive = TRUE)
  saveRDS(ret, "files/6_Neff.RDS")
  cat("Saved: files/6_Neff.RDS\n")
  
  tictoc::toc(quiet = FALSE, func.toc = toc_hours)
  
  return(ret)
}
