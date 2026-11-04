#' @export
#' @describeIn mechGEA Plot migration model results. Consider running with \code{\link{mechGEA}} for the full pipeline.
plot_migration <- function(migration_results, migration_data, patches_xy,
                           mig_plot_threshold = 0.4,
                           mig_plot_topn = 2,
                           mig_plot_binsize = 1){

  patch_from <- patch_to <- p1 <- p2 <- migration_rate <- lon_from <- lon_to <- lat_from <- lat_to <- NULL
  mid_lon <- mid_lat <- mid_lon_bin <-mid_lat_bin <- mig_rate <- lon <- lat <- Distance <- Migration <- NULL

  data <- migration_data; rm(migration_data)
  DistM <- data$DistM
  M <- migration_results; rm(migration_results)

  ## N_ind_per_patch comes from the migration data (already in coordinates order)
  N_ind_per_patch = data$N_table

  ## ---- Build migration network -----------------------------------------------

  diag(M) <- 0
  n = nrow(M)

  migration_lines <- data.frame()
  for (i in 1:n) {
    row_values <- M[i, ]
    row_values[i] <- NA

    topn_indices <- order(row_values, decreasing = TRUE)[1:mig_plot_topn]
    threshold_indices <- which(M[i, ] >= mig_plot_threshold)
    union_indices <- union(topn_indices, threshold_indices)

    for (j in union_indices) {
      migration_lines <- rbind(migration_lines, data.frame(
        patch_from = patches_xy$patch[i],
        patch_to   = patches_xy$patch[j],
        lon_from   = patches_xy$lon[i],
        lat_from   = patches_xy$lat[i],
        lon_to     = patches_xy$lon[j],
        lat_to     = patches_xy$lat[j],
        migration_rate = M[i, j],
        mig_rate = M[i, j]
      ))
    }
  }

  migration_lines <- migration_lines |>
    dplyr::mutate(p1 = pmin(patch_from, patch_to),
                  p2 = pmax(patch_from, patch_to)) |>
    dplyr::distinct(p1, p2, .keep_all = TRUE)

  migration_lines <- migration_lines |> dplyr::arrange(migration_rate)

  ## ---- Reduce dense connections ----------------------------------------------

  migration_lines <- migration_lines |>
    dplyr::mutate(mid_lon = (lon_from + lon_to) / 2,
           mid_lat = (lat_from + lat_to) / 2)

  migration_lines <- migration_lines |>
    dplyr::mutate(mid_lon_bin = round(mid_lon / mig_plot_binsize) * mig_plot_binsize,
           mid_lat_bin = round(mid_lat / mig_plot_binsize) * mig_plot_binsize)

  migration_lines <- migration_lines |>
    dplyr::group_by(mid_lon_bin, mid_lat_bin) |>
    dplyr::filter(migration_rate == max(migration_rate)) |>
    dplyr::ungroup()

  ## ---- Map plot --------------------------------------------------------------

  ## Bounding box: pad by the larger of 20% of the coordinate range or 0.5 degrees,
  ## so that small study areas still show enough geographic context.
  lon_range <- range(patches_xy$lon)
  lat_range <- range(patches_xy$lat)
  pad_lon   <- max(diff(lon_range) * 0.20, 0.5)
  pad_lat   <- max(diff(lat_range) * 0.20, 0.5)
  xlim      <- lon_range + c(-pad_lon, pad_lon)
  ylim      <- lat_range + c(-pad_lat, pad_lat)

  ## Political boundaries (vector data — loads instantly from rnaturalearth)
  world  <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  sf::sf_use_s2(FALSE)   # disable spherical geometry to avoid invalid-edge errors on ne_states
  .suppress_specific_warning(states <- rnaturalearth::ne_states(returnclass = "sf") |>
    sf::st_crop(xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2]), "spatially constant")
  sf::sf_use_s2(TRUE)    # re-enable for everything else

  ## N_ind_per_patch already loaded from data$N_table above

  map <- ggplot2::ggplot() +
    ## Country fill + outlines
    ggplot2::geom_sf(data = world, fill = "grey92", color = "grey40", linewidth = 0.4) +
    ## State / province borders for within-country geographic context
    ggplot2::geom_sf(data = states, fill = NA, color = "grey60", linewidth = 0.25, linetype = "dashed") +
    ## Migration network
    ggplot2::geom_curve(
      data = migration_lines,
      ggplot2::aes(x = lon_from, y = lat_from, xend = lon_to, yend = lat_to,
                   linewidth = mig_rate, color = migration_rate),
      curvature = 0.2, alpha = 0.8,
      show.legend = c(linewidth = FALSE)
    ) +
    ggplot2::scale_linewidth_continuous(trans = "log10", range = c(0.1, 2)) +
    ggplot2::scale_color_gradient(name = "Rel. Migration Rate", low = "yellow", high = "red", trans = "log10") +
    ## Sampling sites
    ggplot2::geom_point(
      data = patches_xy,
      ggplot2::aes(x = lon, y = lat, size = N_ind_per_patch),
      fill = "grey50", color = "black", shape = 21, stroke = 1
    ) +
    ggplot2::scale_size_continuous(range = c(0.5, 4)) +
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::labs(
      x = "Longitude",
      y = "Latitude",
      size = "Individuals \n Sampled"
    ) +
    ggplot2::guides(
      size  = ggplot2::guide_legend(title.position = "top", title.hjust = 0.5),
      color = ggplot2::guide_colorbar(title.position = "top", title.hjust = 0.5),
      linewidth = "none"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    )

  if(!dir.exists("plot")){
    dir.create("plot")
  }
  ggplot2::ggsave("plot/Migration_map.pdf", map, width = 9, height = 8)
  cat("Saved: plot/Migration_map.pdf\n")

  ## ---- Isolation by distance plot --------------------------------------------

  dist_vec <- as.vector(DistM)
  mig_vec <- as.vector(M)
  valid_idx <- which(dist_vec != 0)
  dist_vec <- dist_vec[valid_idx]
  mig_vec <- mig_vec[valid_idx]

  Dist_mig <- ggplot2::ggplot(data.frame(Distance = dist_vec*1e-3, Migration = mig_vec),
                              ggplot2::aes(x = Distance, y = log(Migration))) +
    ggplot2::geom_point(alpha = 0.5, size = 1.2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Distance (km)", y = "log(Rel. Migration Rate)") +
    ggplot2::theme(aspect.ratio = 1) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black")
    )

  ggplot2::ggsave("plot/Isolation_by_distance.pdf", Dist_mig, width = 6, height = 4)
  cat("Saved: plot/Isolation_by_distance.pdf\n")

  ## ---- Combine and save ------------------------------------------------------

  combined_plot <- cowplot::plot_grid(
    Dist_mig, map,
    nrow = 1,
    rel_widths = c(2.8, 3)
  )

  print(combined_plot)

  return(list(IBD = Dist_mig, map = map))
}


#' @export
#' @describeIn mechGEA Plot GEA model results. Consider running with \code{\link{mechGEA}} for the full pipeline.
plot_gea <- function(stan_data,
                     migration_results,
                     beta_draw,
                     beta_summary,
                     rho_kappa_mu,
                     env_names,
                     CI = 0.95,
                     env_colors = NULL){
  Env_var <- Locus <- SNP <- . <- nlogp <- Env <- color <- migration_off <- migration_on <- NULL
  type <- value <- variable <- facet <- muk <- `No migration` <- NULL
  `Mean-field (1st correction)` <- logit <- env_value <- p_logit <- patch <- NULL

  ## ---- Load -----------------------------------------------------------------
  data         <- stan_data; rm(stan_data)
  Migration    <- migration_results; rm(migration_results)

  ## ---- Identify adaptive loci ----------------------------------------------

  result         <- identify_adaptive_loci(beta_draw, beta_summary, data$N_locus,
                                           data$N_env, env_names, CI)
  imp_loci_table <- result$imp_loci_table
  beta_Mlogp     <- result$beta_Mlogp

  if (is.null(imp_loci_table)) {
    inferred <- data.table::data.table(Env_var = character(),
                                       Locus   = integer())
  } else {
    inferred <- imp_loci_table[, .(Env_var, Locus)]
  }

  dir.create("files", showWarnings = FALSE, recursive = TRUE)
  dir.create("plot",  showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(inferred = inferred, CI = CI),
          "files/5_adaptive_loci.RDS")
  cat(sprintf("Inferred %d adaptive (locus, env) pairs at CI = %.2f.\n",
              nrow(inferred), CI))

  ## ---- Env color palette (Okabe-Ito CVD-safe) ------------------------------
  if(is.null(env_colors)){
    okabe_ito <- grDevices::palette.colors(palette = "Okabe-Ito")
    cvd_safe  <- unname(okabe_ito[!names(okabe_ito) %in% c("black", "gray")])
    if (data$N_env > length(cvd_safe))
      cvd_safe <- grDevices::palette.colors(data$N_env, palette = "Polychrome 36",
                                            recycle = FALSE)
    env_colors <- stats::setNames(cvd_safe[seq_len(data$N_env)], env_names)
  } else {
    env_colors <- stats::setNames(env_colors, env_names)
  }


  ## ---- Manhattan plot ------------------------------------------------------
  ymin <- min(beta_Mlogp, na.rm = TRUE)
  ymax <- max(c(beta_Mlogp, -log(1 - CI) + 1), na.rm = TRUE) + 1

  mbMlogp <- data.table::as.data.table(beta_Mlogp)
  colnames(mbMlogp) <- env_names
  mbMlogp[,SNP := 1:.N]
  mbMlogp <- data.table::melt(mbMlogp, id.vars = "SNP", value.name = "nlogp", variable.name = "Env")
  p_man <- ggplot2::ggplot(mbMlogp) +
    ggplot2::geom_point(ggplot2::aes(x = SNP, y = nlogp),
                        color = "grey70", size = .5) +
    ggplot2::xlab("Locus index (candidate set)") +
    ggplot2:: ylab(bquote(-log(p[th])))

  if(nrow(inferred) > 0L) {
    inferred[,Env := match(Env_var, env_names)]
    inferred[,beta := beta_Mlogp[cbind(inferred$Locus, inferred$Env)]]
    inferred[,color := env_colors[inferred$Env_var]]
    inferred

    p_man <- p_man +
      ggplot2::geom_point(data = inferred,
                          ggplot2::aes(x = Locus, y = beta, color = as.character(Env_var)),
                          size = 2, inherit.aes = FALSE) +
      ggplot2::scale_color_manual(breaks = env_names[which(env_names %in% inferred$Env_var)],
                                  values = env_colors[which(env_names %in% inferred$Env_var)]) +
      ggplot2::guides(color = ggplot2::guide_legend("Env"))
    inferred$Env <- inferred$beta <- inferred$color <- NULL
  }

  p_man <- p_man + ggplot2::geom_hline(yintercept = -log(1 - CI),
                                       color = "grey30", inherit.aes = FALSE) +
    ggplot2::theme_bw()

  if(!dir.exists("plot")){
    dir.create("plot")
  }
  ggplot2::ggsave("plot/manhattan.pdf", width = 10, height = 5)
  cat("Saved: plot/manhattan.pdf\n")

  ## ============================================================================
  ## MIGRATION-SELECTION BALANCE for inferred adaptive loci
  ## ============================================================================

  if (nrow(inferred) > 0L) {
    cat("Compiling and exposing migration_off / migration_on ...\n")
    .calc_sel_code <- "
  functions{
    real log_kummer_1F1_pos(real a, real b, real z){
      real log_term = 0;
      real log_sum  = 0;
      for(n in 1:500){
        log_term += log((a + n - 1) * z / ((b + n - 1) * n));
        log_sum   = log_sum_exp(log_sum, log_term);
        if(log_term < log_sum - 30) break;
      }
      return log_sum;
    }
    real log_1F1_full(real a, real b, real z){
      if(z >= 0) return log_kummer_1F1_pos(a, b, z);
      return z + log_kummer_1F1_pos(b - a, b, -z);
    }
    vector pbar0_vec(real mu, real kappa, vector s){
      int N = num_elements(s);
      vector[N] out;
      real a = mu * kappa;
      real b = kappa;
      for(k in 1:N) out[k] = mu * exp(log_1F1_full(a + 1, b + 1, s[k])
                                    - log_1F1_full(a,     b,     s[k]));
      return out;
    }
    matrix migration_off(int N_patch, int N_iter,
                         vector kappa, vector rho, vector mu,
                         matrix migration, matrix s_coef){
      matrix [N_patch, N_iter] p_off;
      for(i in 1:N_iter) p_off[,i] = pbar0_vec(mu[i], kappa[i], s_coef[,i]);
      return p_off;
    }
    matrix migration_on(int N_patch, int N_iter,
                        vector kappa, vector rho, vector mu,
                        matrix migration, matrix s_coef){
      matrix [N_patch, N_iter] p_on;
      vector[N_patch] diag_M = diagonal(migration);
      vector[N_patch] S      = -diag_M;
      for(i in 1:N_iter){
        vector[N_patch] kappa_star      = kappa[i] + rho[i] * S;
        vector[N_patch] pbar0           = pbar0_vec(mu[i], kappa[i], s_coef[,i]);
        vector[N_patch] M_offdiag_pbar0 = migration * pbar0 - diag_M .* pbar0;
        vector[N_patch] mu_star = (mu[i] * kappa[i]
                                    + rho[i] * M_offdiag_pbar0) ./ kappa_star;
        for(k in 1:N_patch){
          real ak = mu_star[k] * kappa_star[k];
          real bk = kappa_star[k];
          p_on[k,i] = mu_star[k] * exp(log_1F1_full(ak + 1, bk + 1, s_coef[k,i])
                                     - log_1F1_full(ak,     bk,     s_coef[k,i]));
        }
      }
      return p_on;
    }
  }
  "
    rstan::expose_stan_functions(rstan::stanc(model_code = .calc_sel_code))

    N_locus_all <- data$N_locus
    N_iter      <- nrow(beta_draw)
    .rho   <- rho_kappa_mu[, 1:N_locus_all]
    .kappa <- rho_kappa_mu[, (N_locus_all + 1):(2 * N_locus_all)]
    .mu    <- rho_kappa_mu[, (2 * N_locus_all + 1):(3 * N_locus_all)]

    bal_table     <- as.data.frame(inferred[order(Locus)])
    bal_table$Env <- match(bal_table$Env_var, env_names)

    p_mig_off <- array(0, c(data$N_patch, N_iter, nrow(bal_table)))
    p_mig_on  <- array(0, c(data$N_patch, N_iter, nrow(bal_table)))

    cat("Computing per-deme means for each adaptive locus ...\n")
    for (i in seq_len(nrow(bal_table))) {
      l <- bal_table$Locus[i]
      e <- bal_table$Env[i]
      beta_col <- (l - 1) * data$N_env + e
      selection <- matrix(beta_draw[, beta_col],
                          nrow = data$N_patch, ncol = N_iter, byrow = TRUE) *
        matrix(as.vector(unlist(data$Env[, e])),
               nrow = data$N_patch, ncol = N_iter)
      p_mig_off[, , i] <- migration_off(data$N_patch, N_iter,
                                        .kappa[, l], .rho[, l], .mu[, l],
                                        Migration, selection)
      p_mig_on [, , i] <- migration_on (data$N_patch, N_iter,
                                        .kappa[, l], .rho[, l], .mu[, l],
                                        Migration, selection)
      if (i %% 5L == 0L || i == nrow(bal_table))
        cat(sprintf("  locus %2d / %2d done\n", i, nrow(bal_table)))
    }

    # prepare plot data

    col_off    <- "#0072B2"
    col_on     <- "#E69F00"
    col_arrow1 <- "grey25"


    melt_ptable <- function(x){
      x <- data.table::as.data.table(x)
      colnames(x) <- paste0("Locus ", bal_table$Locus)
      x <- cbind(x, patch = 1:nrow(data$Env))
      x <- data.table::melt(data.table::as.data.table(x), id.vars = "patch")
    }

    env_mat <- melt_ptable(data$Env[, bal_table$Env])
    env_mat[,Env := rep(bal_table$Env, each = nrow(data$Env))]
    env_mat[,Env := env_names[Env]]
    y_off <- melt_ptable(apply(p_mig_off, 3, rowMeans))
    y_on <- melt_ptable(apply(p_mig_on, 3, rowMeans))

    y_off[,type := "No migration"]
    y_on[,type := "Mean-field  (1st correction)"]

    y_off$env_value <- env_mat$value
    y_on$env_value <- env_mat$value
    y_off$Env <- env_mat$Env
    y_on$Env <- env_mat$Env

    msb_data <- rbind(y_on, y_off)
    msb_data <- dplyr::rename(msb_data,
                              muk = value,
                              Locus = variable)

    cast_msb <- data.table::dcast(msb_data, patch + Locus + env_value + Env ~ type, value.var = "muk")


    # plot a multi-page pdf and generate per-page ggplots.
    panels <- unique(msb_data[,.(Locus, Env)])
    panels <- dplyr::arrange(panels, Locus, Env)

    panels_per_page <- 16
    total_pages <- ceiling(nrow(panels)/panels_per_page)
    data.table::set(panels, i = 1:nrow(panels), j = "page", value = rep(1:total_pages, each = panels_per_page)[1:nrow(panels)])
    panels[,facet := paste0(Locus, "_", Env)]

    msb_plots <- vector("list", length(unique(panels$page)))
    names(msb_plots) <- paste0("page_", unique(panels$page))

    grDevices::pdf("plot/migration_sel_balance.pdf", width = 8.5, height = 9.4)

    for(i in 1:length(msb_plots)){
      tpd <- msb_data[paste0(Locus, "_", Env) %in% panels[panels$page == i,]$facet]
      tcast <- cast_msb[paste0(Locus, "_", Env) %in% panels[panels$page == i,]$facet]

      msb_plots[[i]] <- ggplot2::ggplot(tpd, ggplot2::aes(x = env_value, y = muk, color = type)) +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_segment(data = tcast, ggplot2::aes(y = `No migration`,
                                                         yend = `Mean-field  (1st correction)`,
                                                         x = env_value, xend = env_value),
                              arrow = grid::arrow(length = grid::unit(6, "points")),
                              inherit.aes = FALSE) +
        ggh4x::facet_wrap2(Locus ~ Env, scales = "free_x", strip = ggh4x::strip_split(position = c("top", "bottom"))) +
        ggplot2::scale_color_manual(breaks = c("No migration", "Mean-field  (1st correction)"),
                                    values = c(col_off, col_on)) +
        ggplot2::theme_bw() +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::theme(legend.title = ggplot2::element_blank(),
                       strip.background = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank(),
                       strip.placement = "outside") +
        ggplot2::ylab(bquote('Per-deme mean' ~ mu[k])) +
        ggplot2::ggtitle("Migration-Selection balance for adaptive loci")

      print(msb_plots[[i]])
    }
    grDevices::dev.off()
    cat("Saved: plot/migration_sel_balance.pdf\n")


  } else {
    cat("No inferred adaptive loci -- skipping migration-sel balance plot.\n")
    msb_plots <- NULL
  }

  ## ============================================================================
  ## GE-ASSOCIATION via posterior_check.stan
  ## ============================================================================

  if (nrow(inferred) == 0L) {
    cat("No inferred adaptive loci -- skipping posterior_check.stan.\n")
    asso_plots <- NULL
  } else {
    inferred      <- inferred[order(Locus)]
    imp_env_index <- match(inferred$Env_var, env_names)

    imp_loci_data <- list(
      N_patch  = data$N_patch,
      N_locus  = nrow(inferred),
      G_eff    = data$G_call[inferred$Locus, , drop = FALSE],
      TwoN_eff = data$TwoN  [inferred$Locus, , drop = FALSE],
      Env      = matrix(NA_real_, nrow = nrow(inferred), ncol = data$N_patch)
    )
    for (i in seq_len(nrow(inferred))) {
      imp_loci_data$Env[i, ] <- data$Env[, imp_env_index[i]]
    }

    STAN_MODEL <- .source_model("posterior_check.stan")
    if(!dir.exists("./stan")){dir.create("stan")}
    if(!file.exists("stan/posterior_check.stan")){
      file.copy(STAN_MODEL, "./stan/posterior_check.stan")
    }


    mod_pc <- cmdstanr::cmdstan_model(
      "stan/posterior_check.stan",
      cpp_options = list(stan_threads = TRUE)
    )

    # DEV NOTE: is it worth allowing these to change?
    n_cores_pc <- parallel::detectCores()
    chain_pc   <- 2
    threads_pc <- max(1L, ceiling(n_cores_pc / chain_pc))

    fit_pc <- mod_pc$sample(
      data              = imp_loci_data,
      chains            = chain_pc,
      seed              = 123,
      parallel_chains   = chain_pc,
      threads_per_chain = threads_pc,
      refresh           = 1000,
      iter_warmup       = 350,
      iter_sampling     = 800,
      show_messages     = TRUE
    )

    p_mean <- matrix(fit_pc$summary(variables = "p")$mean,
                     nrow = imp_loci_data$N_locus, ncol = imp_loci_data$N_patch)

    # format data for plotting
    melt_pmean <- function(x){
      x <- data.table::as.data.table(t(x))
      colnames(x) <- paste0("Locus ", inferred$Locus)
      x[,patch := paste0("patch_", 1:nrow(data$Env))]
      return(data.table::melt(x, id.vars = "patch"))
    }

    p_mean <- melt_pmean(p_mean)
    env <- melt_pmean(imp_loci_data$Env)
    env[,Env := rep(inferred$Env_var, each = nrow(data$Env))]
    p_mean[,logit := LaplacesDemon::logit(value)]


    p_mean$env_value <- env$value
    p_mean$Env <- env$Env

    ymin_l <- min(p_mean$logit)
    ymax_l <- max(p_mean$logit)

    p_mean <- dplyr::rename(p_mean, p_logit = logit, Locus = variable)

    panels <- unique(p_mean[,.(Locus, Env)])
    panels <- dplyr::arrange(panels, Locus, Env)


    # plot the associations
    panels_per_page <- 16
    total_pages <- ceiling(nrow(panels)/panels_per_page)
    data.table::set(panels, i = 1:nrow(panels), j = "page", value = rep(1:total_pages, each = panels_per_page)[1:nrow(panels)])
    panels[,facet := paste0(Locus, "_", Env)]

    asso_plots <- vector("list", length(unique(panels$page)))
    names(asso_plots) <- paste0("page_", unique(panels$page))

    grDevices::pdf("plot/GE_association.pdf", width = 8, height = 9.4)

    for(i in 1:length(asso_plots)){
      tpd <- p_mean[paste0(Locus, "_", Env) %in% panels[panels$page == i,]$facet]

      asso_plots[[i]] <- ggplot2::ggplot(tpd, ggplot2::aes(x = env_value, y = p_logit, color = Env)) +
        ggplot2::geom_point() +
        ggplot2::geom_smooth(method = "lm", se = FALSE) +
        ggh4x::facet_wrap2(Locus ~ Env, scales = "free_x", strip = ggh4x::strip_split(position = c("top", "bottom"))) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       strip.background = ggplot2::element_blank(),
                       strip.placement = "outside",
                       legend.position = "none") +
        ggplot2::scale_y_continuous(limits = c(ymin_l, ymax_l)) +
        ggplot2::scale_color_manual(breaks = unique(as.character(inferred$Env_var)),
                                    values = unique(env_colors[inferred$Env_var])) +
        ggplot2::ylab(paste0("logit(p)"))

      .suppress_specific_warning(print(asso_plots[[i]]), "outside the scale range")

    }
    grDevices::dev.off()
    cat("Saved: plot/GE_association.pdf\n")
  }

  return(list(SNP_pvalues = mbMlogp,
              manhattan = p_man,
              migration_selection_balance = msb_plots,
              association_plots = asso_plots))
}

#' Generate a manhattan plot one or more mechGEA outputs.
#'
#' Creates a ggplot-based manhattan plot, where chromosomes/scaffolds/etc are
#' concatenated along the x-axis. Works with \code{\link{mechGEA}} SNP_pvalue
#' outputs; take either a single output or a list of outputs from multiple
#' chunks (see details).
#'
#' @param SNP_pvalues data.frame or similar like that produced by
#'   \code{\link{mechGEA}} in \code{$SNP_pvalues} OR the output from
#'   \code{\link{mechGEA}} \ \code{\link{run_stan_GEA_model}} OR a list of a
#'   either of the two (mixing OK). If the former, contains three columns named
#'   "SNP", "Env", and "nlogp" containing SNP IDs, names for environmental
#'   effects, and -log p-values for associations.
#' @param SNP_metadata data.frame or similar containing metadata on SNPs OR a
#'   path to a ".beagle" or ".beagle.gz" file containing SNP metadata OR a list
#'   of a either of the two (mixing OK). If the former, must have three columns
#'   named "SNP", "CHROM", and "POS" containing SNP IDs, chromosomes, and
#'   positions for each SNP. Each SNP ID in the \code{SNP_pvalues} "SNP" column
#'   must be represented in the "SNP" column of \code{SNP_metadata}. Assumes
#'   sorting is in the same order as in \code{SNP_pvalues}, across lists or not.
#' @param vlines character (color) or FALSE, default FALSE. If a color, vertical
#'   separator lines will be drawn between each chromosome. Widths controlled by
#'   \code{vline_width}.
#' @param vline_width numeric, default 2. Width of chromosome separator lines.
#'   Ignored if \code{vlines} is FALSE.
#' @param significant numeric, default NULL. Value at which a line will be drawn
#'   designating significant SNPs. If highlight = "significant", SNPs above this
#'   level will also be labeled.
#' @param suggestive numeric, default NULL. Value at which a line will be drawn
#'   designating suggestive SNPs. If highlight = "suggestive", SNPs above this
#'   level will also be labeled.
#' @param colors character, default c("black", "slategray3"). Colors to
#'   alternate across chromosomes.
#' @param chr_order character, default NULL. If provided, an ordered vector of
#'   chromosome/scaffold/etc names by which to sort output.
#' @param abbreviate_labels numeric or FALSE, default FALSE. If a numeric value,
#'   x-axis chromosome names will be abbreviated using
#'   \code{\link[base]{abbreviate}}, with each abbreviated label having the
#'   minimum length specified. Helpful when chromosome/scaffold/etc names are
#'   very long.
#'
#' @author William Hemstrom
#' @export
#' @examples
#' \dontrun{
#' # run a model
#' beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
#' phred <- beagle2pl(beagle, lines = 1:150)
#' input_data <- prepare_metadata(sample_metadata)
#' res <- mechGEA(phred,
#'                patch_info = input_data$patch_info,
#'                patches_xy = input_data$patches_xy,
#'                migration_loci = 50,
#'                env = input_data$env,
#'                iter_sampling_gea = 100,
#'                iter_warmup_gea = 200,
#'                iter_sampling_mig = 100,
#'                iter_warmup_mig = 200)
#'
#' # plot results
#' plot_manhattan(res, beagle)
#' ## or
#' plot_manhattan(res$SNP_pvalues, beagle)
#' # or could provide a data.frame of SNP metadata
#'
#'
#' # combine across result runs. Note that using chunk_beagle() is easier, than
#' # manually splitting, but for the example we're just making two chunks
#' beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
#' system(paste0("zcat ", beagle2,
#'               " | head -n ", 301, " | gzip > test_sub.beagle.gz"))
#' b <- fread(" test_sub.beagle.gz")
#' b <- b[151:300,] # the second batch of lines
#' fwrite(b, "test_sub.beagle.gz", sep = "\t",
#'        col.names = TRUE, row.names = FALSE)
#'
#' phred <- beagle2pl("test_sub.beagle.gz")
#' res2 <- mechGEA(phred,
#'                 patch_info = input_data$patch_info,
#'                 patches_xy = input_data$patches_xy,
#'                 migration_loci = 50,
#'                 env = input_data$env,
#'                 iter_sampling_gea = 100,
#'                 iter_warmup_gea = 200,
#'                 iter_sampling_mig = 100,
#'                 iter_warmup_mig = 200)
#'
#' # plot results
#' plot_manhattan(list(res, res2), list(beagle, "test_sub.beagle.gz"))
#' # Any combination of res/res$SNP_pvalues and filepaths/tables of metadata
#' # are allowed in lists. Assumes everything IS IN THE SAME ORDER!
#' }
#'
plot_manhattan <- function(SNP_pvalues,
                           SNP_metadata,
                           vlines = FALSE,
                           vline_width = .25,
                           significant = 5, suggestive = 3,
                           colors = c("black", "slategray3"),
                           chr_order = NULL,
                           abbreviate_labels = FALSE,
                           env_colors = NULL){

  cum.bp <- cum.start <- cum.end <- y <- start <- end <- med <- POS <- IDfac <- SNP <-  NULL
  #=============process input data==============================
  # processing the input p-values
  num_pvals <- 0
  if(is.list(SNP_pvalues) & !is.data.frame(SNP_pvalues)){
    if("SNP_pvalues" %in% names(SNP_pvalues)){
      SNP_pvalues <- SNP_pvalues$SNP_pvalues
    } else {
      combined_pvalues <- vector("list", length(SNP_pvalues))
      for(i in 1:length(SNP_pvalues)){

        if(!is.data.frame(SNP_pvalues[[i]])){
          if("SNP_pvalues" %notin% names(SNP_pvalues[[i]])){
            stop(paste0("Could not find data in SNP_pvalues[[", i, "]]. Is this data either a list of data.frames/tables with SNP_pvalues or a list of outputs from mechGEA?"))
          }
          combined_pvalues[[i]] <- SNP_pvalues[[i]][["SNP_pvalues"]]
        } else {
          combined_pvalues[[i]] <- SNP_pvalues[[i]]
        }

        combined_pvalues[[i]][["SNP"]] <- combined_pvalues[[i]][["SNP"]] + num_pvals
        num_pvals <- num_pvals + length(unique(combined_pvalues[[i]][["SNP"]]))
        dreamerr::check_set_arg(combined_pvalues[[i]], .type = "data.frame ncol(3)", .message = paste0("Object in either SNP_pvalues[[", i, "]] or SNP_pvalues[[", i, "]][['SNP_pvalues']] does not appear to be a three column data.frame/table.\n"))
      }
    }
  }

  # if a beagle/list of beagles is instead handed over for the metadata, process those for CHR/position.
  ## beagle meta read subfunc
  beagle_read_meta <- function(file){
    if(file.exists(file)) {
      if(grepl("\\.beagle\\.gz$", file) | grepl("\\.beagle$", file)){
        res <- data.table::fread(file, select = 1)
        tcol <- colnames(res)[1]
        res[, c("CHROM", "POS") := data.table::tstrsplit(get(tcol), "_")]
        res[,POS := as.numeric(POS)]
        res[[1]] <- 1:nrow(res)
      } else {stop(paste0("File: ", file, " does not end in '.beagle' or '.beagle.gz'.\n"))}
    } else {stop(paste0("Could not locate file: ", file, "\n"))}
    colnames(res) <- c("SNP", "CHROM", "POS")
    return(res)
  }


  ## read it in
  num_meta <- 0
  if(is.list(SNP_metadata) & !is.data.frame(SNP_metadata)){
    combined_meta <- vector("list", length(SNP_metadata))

    for(i in 1:length(SNP_metadata)){
      if(!is.data.frame(SNP_metadata[[i]])){
        combined_meta[[i]] <- beagle_read_meta(SNP_metadata[[i]])
      } else {
        combined_meta[[i]] <- SNP_metadata[[i]]
      }

      combined_meta[[i]][["SNP"]] <- num_meta + combined_meta[[i]][["SNP"]]
      num_meta <- num_meta + nrow(combined_meta[[i]])
    }
  } else if(!is.data.frame(SNP_metadata)){
    SNP_metadata <- beagle_read_meta(SNP_metadata)
  }

  # merge if needed
  if(exists("combined_pvalues")){SNP_pvalues <- data.table::rbindlist(combined_pvalues, fill = TRUE)}
  if(exists("combined_meta")){SNP_metadata <- rbindlist(combined_meta, fill = TRUE)}

  #=============sanity checks==============================
  dreamerr::check_set_arg(SNP_pvalues, SNP_metadata, .type = "data.frame ncol(3)")
  SNP_pvalues <- data.table::as.data.table(SNP_pvalues)
  SNP_metadata <- data.table::as.data.table(SNP_metadata)
  msg <- character()


  # columns OK?
  expected_cols <- c("SNP", "CHROM", "POS")
  missing_cols <- expected_cols %notin% colnames(SNP_metadata)
  if(any(missing_cols)){
    msg <- c(msg, paste0("Some expected column names not found in SNP metadata. Missing: \n\t", paste0(expected_cols[missing_cols], collapse = "\n\t")))
  }

  expected_cols <- c("SNP", "Env", "nlogp")
  missing_cols <- expected_cols %notin% colnames(SNP_pvalues)
  if(any(missing_cols)){
    msg <- c(msg, paste0("Some expected column names not found in SNP pvalues. Missing: \n\t", paste0(expected_cols[missing_cols], collapse = "\n\t")))
  }

  # all SNPs accounted for?
  bad_snps <- which(SNP_pvalues$SNP %notin% SNP_metadata$SNP)
  if(length(bad_snps) > 0){
    msg <- c(msg, paste0("Some SNP IDs in pvals not found in metadata. Number not found: ", length(bad_snps), "\n"))
  }

  if(length(msg) > 0){
    stop(msg, collapse = "\n")
  }

  SNP_pvalues <- merge(SNP_pvalues, SNP_metadata)

  stats <- SNP_pvalues; rm(SNP_pvalues)

  plot_var <- "nlogp"; chr <- "CHROM"; bp <- "POS"

  if(nrow(stats) == 0){
    stop("No non-NA data found.\n")
  }


  #=============figure out adjusted positions on x axis================
  # fetch chromosome info
  if(is.factor(stats[[chr]])){
    stats[[chr]] <- as.character(stats[[chr]])
  }

  if(!is.null(chr_order)){
    stats[[chr]] <- factor(stats[[chr]], levels = chr_order)
  }

  # get adjusted x axis positions and figure out where chromsome tick marks should be placed.
  chr.info <- tapply(stats[[bp]], stats[[chr]], max) # chromosome lengths
  chr.centers <- (chr.info - tapply(stats[[bp]], stats[[chr]], min))/2 # chromosome centers
  chr.info <- cumsum(chr.info) # cumulative lengths
  cum.bp <- c(0, chr.info[-length(chr.info)]) # correct cummulative lengths
  names(cum.bp) <- names(chr.info) # rename
  cum.chr.centers <- cum.bp + chr.centers # cumulative chromosome centers
  stats$start <- cum.bp[match(stats[[chr]], names(cum.bp))]
  stats[,cum.bp := start + POS]

  #============produce the plot========
  if(is.null(env_colors)){
    env_colors <- pal_picker(length(unique(stats[nlogp >= significant,]$Env)))
  }

  p <- ggplot2::ggplot(stats, ggplot2::aes(x = cum.bp, y = nlogp, color = as.factor(CHROM))) +
    ggplot2::geom_point(show.legend = FALSE, size = 1) +
    ggplot2::scale_color_manual(values = rep(c(colors), length(cum.chr.centers))) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(data = stats[nlogp >= significant,], ggplot2::aes(x = cum.bp, y = nlogp, color = as.factor(Env)),
                        size = 3) +
    ggplot2::scale_color_manual(values = env_colors)

  p <- p +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(hjust = 0.01)) +
    ggplot2::xlab("Chromosome")

  p <- p + ggplot2::ylab(bquote(-log[10] ~ italic("p")))
  p <- p + ggplot2::guides(color = ggplot2::guide_legend("Env"))


  #=============adjust the plot========
  # significant/suggestive lines
  if(!is.null(significant)){
    p <- p + ggplot2::geom_hline(yintercept = significant, color = "#DDAA33")
  }
  if(!is.null(suggestive)){
    p <- p + ggplot2::geom_hline(yintercept = suggestive, color = "#004488")
  }


  # abbreviate
  if(!isFALSE(abbreviate_labels)){
    p <- p + ggplot2::scale_x_continuous(label = abbreviate(names(cum.chr.centers),
                                                            minlength = abbreviate_labels),
                                         breaks = cum.chr.centers, minor_breaks = NULL)
  }
  else{
    p <- p + ggplot2::scale_x_continuous(label = names(cum.chr.centers),
                                         breaks = cum.chr.centers, minor_breaks = NULL)
  }


  if(!isFALSE(vlines)){

    cum.start <- cum.chr.centers + tapply(stats[[bp]], stats[[chr]], function(x) max(x)/2)
    p <- p + ggplot2::geom_vline(xintercept = as.numeric(cum.start), color = vlines, linewidth = vline_width)
  }

 return(p)
}
