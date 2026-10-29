#' @export
#' @describeIn mechGEA Plot migration model results. Consider running with \code{\link{mechGEA}} for the full pipeline.
plot_migration <- function(migration_results, migration_data, patches_xy,
                           mig_plot_threshold = 0.4,
                           mig_plot_topn = 2,
                           mig_plot_binsize = 1){

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

  # DEV NOTES: Plots should probably all be swapped over to ggplots and the raw data
  # returned to allow for much easier customization of plots!

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
    env_colors <- setNames(cvd_safe[seq_len(data$N_env)], env_names)
  } else {
    env_colors <- setNames(env_colors, env_names)
  }


  ## ---- Manhattan plot ------------------------------------------------------
  ymin <- min(beta_Mlogp, na.rm = TRUE)
  ymax <- max(c(beta_Mlogp, -log(1 - CI) + 1), na.rm = TRUE) + 1

  mbMlogp <- as.data.table(beta_Mlogp)
  colnames(mbMlogp) <- env_names
  mbMlogp[,SNP := 1:.N]
  mbMlogp <- melt(mbMlogp, id.vars = "SNP", value.name = "nlogp", variable.name = "Env")
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
      x <- as.data.table(x)
      colnames(x) <- paste0("Locus ", bal_table$Locus)
      x <- cbind(x, patch = 1:nrow(data$Env))
      x <- melt(as.data.table(x), id.vars = "patch")
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

    cast_msb <- dcast(msb_data, patch + Locus + env_value + Env ~ type, value.var = "muk")


    # plot a multi-page pdf and generate per-page ggplots.
    panels <- unique(msb_data[,.(Locus, Env)])
    panels <- dplyr::arrange(panels, Locus, Env)

    panels_per_page <- 16
    total_pages <- ceiling(nrow(panels)/panels_per_page)
    set(panels, i = 1:nrow(panels), j = "page", value = rep(1:total_pages, each = panels_per_page)[1:nrow(panels)])
    panels[,facet := paste0(Locus, "_", Env)]

    msb_plots <- vector("list", length(unique(panels$page)))
    names(msb_plots) <- paste0("page_", unique(panels$page))

    pdf("plot/migration_sel_balance.pdf", width = 8.5, height = 9.4)

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
    dev.off()
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


    mod_pc <- cmdstan_model(
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
      x <- as.data.table(t(x))
      colnames(x) <- paste0("Locus ", inferred$Locus)
      x[,patch := paste0("patch_", 1:nrow(data$Env))]
      return(melt(x, id.vars = "patch"))
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
    set(panels, i = 1:nrow(panels), j = "page", value = rep(1:total_pages, each = panels_per_page)[1:nrow(panels)])
    panels[,facet := paste0(Locus, "_", Env)]

    asso_plots <- vector("list", length(unique(panels$page)))
    names(asso_plots) <- paste0("page_", unique(panels$page))

    pdf("plot/GE_association.pdf", width = 8, height = 9.4)

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
    dev.off()
    cat("Saved: plot/GE_association.pdf\n")
  }

  return(list(SNP_pvalues = mbMlogp,
              manhattan = p_man,
              migration_selection_balance = msb_plots,
              association_plots = asso_plots))
}
