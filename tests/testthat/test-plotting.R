test_that("gea plots",{
  skip_on_cran(); skip_on_ci()

  beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
  phred <- beagle2pl(beagle, lines = 1:150)

  input_data <- prepare_metadata(sample_metadata)

  # running these tests are tricky, because they require dependencies that
  # are tough to install without the env script. As a result, we will run this
  # on an interactive session and save results. Tests will only check pre-run
  # results. This DOES mean we need to manually run the tests when updates can
  # cause issues, which stinks.
  tfile <- test_path("testdata", "small_example_results.RDS")
  if(interactive()){
    set.seed(42)
    res <- mechGEA(phred,
                   patch_info = input_data$patch_info,
                   patches_xy = input_data$patches_xy,
                   migration_loci = 50,
                   env = input_data$env,
                   iter_sampling_gea = 100,
                   iter_warmup_gea = 200,
                   iter_sampling_mig = 100,
                   iter_warmup_mig = 200)
    saveRDS(res, tfile)
  } else {
    skip_if_not(file.exists(tfile),
                "test data not available")

    res <- readRDS(tfile)
  }

  p <- plot_gea(stan_data = res$stan_and_migration_data$stan_data,
                migration_results = res$migration_results,
                beta_draw = res$GEA_results$beta_draw,
                beta_summary = res$GEA_results$beta_summary,
                rho_kappa_mu =  res$GEA_results$rho_kappa_mu,
                env_names = colnames(input_data$env)[-1],
                CI = 0.95)

  # are there reasonable plots?
  expect_true(methods::is(p$manhattan, "ggplot"))
  expect_true(all(unlist(lapply(p$migration_selection_balance, methods::is, class = "ggplot"))))
  expect_true(all(unlist(lapply(p$association_plots, methods::is, class = "ggplot"))))

  expect_equal(sum(p$manhattan@data$nlogp >= 4), 5)
  expect_equal(nrow(ggplot2::ggplot_build(p$migration_selection_balance[[1]])$layout$layout),
               nrow(ggplot2::ggplot_build(p$association_plots[[1]])$layout$layout))

  # pvals
  pvals <- p$SNP_pvalues
  senvgrid <- expand.grid(Locus = 1:150, env = names(input_data$env)[-1])
  senvgrid <- data.table::as.data.table(senvgrid)
  senvgrid[,facet := paste0(Locus, "_", env)]
  expect_identical(sort(paste0(pvals$SNP, "_", pvals$Env)),
                   sort(senvgrid$facet)) # everything accounted for?
  expect_equal(sum(pvals$nlogp >= 4), 5)


})

test_that("migration_plots",{

  beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
  phred <- beagle2pl(beagle, lines = 1:150)

  input_data <- prepare_metadata(sample_metadata)

  # running these tests are tricky, because they require dependencies that
  # are tough to install without the env script. As a result, we will run this
  # on an interactive session and save results. Tests will only check pre-run
  # results. This DOES mean we need to manually run the tests when updates can
  # cause issues, which stinks.
  tfile <- test_path("testdata", "small_example_results.RDS")
  if(interactive()){
    res <- mechGEA(phred,
                   patch_info = input_data$patch_info,
                   patches_xy = input_data$patches_xy,
                   migration_loci = 50,
                   env = input_data$env,
                   iter_sampling_mig = 100,
                   iter_warmup_mig = 200,
                   stop_after_migration = TRUE)
  } else {
    skip_if_not(file.exists(tfile),
                "test data not available")

    res <- readRDS(tfile)
  }


  p <- plot_migration(migration_results = res$migration_results,
                 migration_data = res$stan_and_migration_data$migration_data,
                 patches_xy = input_data$patches_xy,
                 mig_plot_threshold = 0.4,
                 mig_plot_topn = 2,
                 mig_plot_binsize = 1)

  expect_true(methods::is(p$IBD, "ggplot"))
  expect_true(methods::is(p$map, "ggplot"))
})

test_that("manhattan_plots", {
  # with one outfile; multiple outfiles will simply merge them.
  tfile <- test_path("testdata", "small_example_results.RDS")
  res <- readRDS(tfile)


  # position info, which we generate here but is usally provided
  snp_metadata <- data.table::data.table(SNP = unique(res$SNP_pvalues$SNP),
                                         CHROM = sort(rep(paste0("CHROM.", 1:10), length.out = length(unique(res$SNP_pvalues$SNP)))))
  snp_metadata[, POS := 1:.N, by = CHROM]

  man <- plot_manhattan(res$SNP_pvalues, SNP_metadata = snp_metadata)
  man2 <- plot_manhattan(res, SNP_metadata = snp_metadata)

  expect_true(ggplot2::is_ggplot(man)) # other tests can come later for basic plotting
  expect_identical(man, man2)

  ## test metadata matching from beagle
  beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
  system(paste0("zcat ", beagle, " | head -n ", 151, " | gzip > test_sub.beagle.gz"))
  bf <- fread("test_sub.beagle.gz")
  bf$marker <- paste0(snp_metadata$CHROM, "_", snp_metadata$POS)
  fwrite(bf, "test_sub.beagle.gz", sep = "\t", col.names = TRUE, row.names = FALSE)
  man3 <- plot_manhattan(res$SNP_pvalues, "test_sub.beagle.gz", chr_order = paste0("CHROM.", 1:10))
  expect_true(ggplot2::is_ggplot(man3)) # test the data in a minute
  file.remove("test_sub.beagle.gz")

  ## test that reordering chrs works
  man1_reorder <- plot_manhattan(res$SNP_pvalues, SNP_metadata = snp_metadata, chr_order = paste0("CHROM.", 1:10))
  xlabs <- ggplot2::ggplot_build(man1_reorder)
  expect_equal(xlabs$layout$panel_params[[1]]$x$get_labels(), paste0("CHROM.", 1:10))
  ### test that the values are sorted the same as when we did metadata from beagle, above
  expect_equal(dplyr::arrange(man3@data[,.(CHROM, POS, start, cum.bp, nlogp, Env)], CHROM, POS, Env),
               dplyr::arrange(man1_reorder@data[,.(CHROM, POS, start, cum.bp, nlogp, Env)], CHROM, POS, Env)) # same data means same plot.

  #======================multiple plots -- there's a TON to check on metadata matching.==================
  skip_on_cran(); skip_on_ci()

  # prep--run another set
  beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
  phred <- beagle2pl(beagle, lines = 151:300)

  input_data <- prepare_metadata(sample_metadata)

  res2 <- mechGEA(phred,
                 patch_info = input_data$patch_info,
                 patches_xy = input_data$patches_xy,
                 migration_loci = 50,
                 env = input_data$env,
                 iter_sampling_gea = 100,
                 iter_warmup_gea = 200,
                 iter_sampling_mig = 100,
                 iter_warmup_mig = 200)

  snp_metadata2 <- data.table::data.table(SNP = unique(res2$SNP_pvalues$SNP),
                                          CHROM = sort(rep(paste0("CHROM.", 1:10), length.out = length(unique(res$SNP_pvalues$SNP)))))
  snp_metadata2[, POS := (1:.N) + 15, by = CHROM]

  # basic plot, direct data for each.
  man <- plot_manhattan(list(res$SNP_pvalues, res2$SNP_pvalues), SNP_metadata = list(snp_metadata,
                                                                                     snp_metadata2))
  expect_equal(length(unique(paste0(man@data$CHROM, "_", man@data$POS))), 300)
  expect_true(all(man@data$POS[man@data$SNP >= 151] >= 16))

  # mix for pvals
  man2 <- plot_manhattan(list(res$SNP_pvalues, res2), SNP_metadata = list(snp_metadata,
                                                                          snp_metadata2))
  expect_equal(man@data, man2@data)

  # mix for metadata sources
  beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
  system(paste0("zcat ", beagle, " | head -n ", 301, " | gzip > test_sub.beagle.gz"))
  bf <- fread("test_sub.beagle.gz")
  bf <- bf[151:300,]
  bf$marker <- paste0(snp_metadata2$CHROM, "_", snp_metadata2$POS)
  fwrite(bf, "test_sub.beagle.gz", sep = "\t", col.names = TRUE, row.names = FALSE)

  man2 <- plot_manhattan(list(res, res2$SNP_pvalues), SNP_metadata = list(snp_metadata,
                                                              "test_sub.beagle.gz"))
  expect_equal(dplyr::arrange(man@data[,.(CHROM, POS, start, cum.bp, nlogp, Env)], CHROM, POS, Env),
               dplyr::arrange(man2@data[,.(CHROM, POS, start, cum.bp, nlogp, Env)], CHROM, POS, Env)) # same data means same plot.


  # a combined data.frame for only one
  r2t <-  copy(res2$SNP_pvalues)
  r2t[,SNP := SNP + 150]
  man3 <- plot_manhattan(rbind(res$SNP_pvalues, r2t), SNP_metadata = list(snp_metadata, "test_sub.beagle.gz"))
  expect_equal(dplyr::arrange(man2@data[,.(CHROM, POS, start, cum.bp, nlogp, Env)], CHROM, POS, Env),
               dplyr::arrange(man3@data[,.(CHROM, POS, start, cum.bp, nlogp, Env)], CHROM, POS, Env)) # same data means same plot.


  # a split SNP data but a combined table of metadata
  r2t <- copy(snp_metadata2)
  r2t[,SNP := SNP + 150]
  man3 <- plot_manhattan(list(res, res2$SNP_pvalues), SNP_metadata = rbind(snp_metadata, r2t))
  expect_equal(dplyr::arrange(man2@data[,.(CHROM, POS, start, cum.bp, nlogp, Env)], CHROM, POS, Env),
               dplyr::arrange(man3@data[,.(CHROM, POS, start, cum.bp, nlogp, Env)], CHROM, POS, Env)) # same data means same plot.

  # lastly, split SNP data but a combined BEAGLE file
  system(paste0("zcat ", beagle, " | head -n ", 301, " | gzip > test_sub.beagle.gz"))
  bf <- fread("test_sub.beagle.gz")
  bf$marker <- c(paste0(snp_metadata$CHROM, "_", snp_metadata$POS), paste0(snp_metadata2$CHROM, "_", snp_metadata2$POS))
  fwrite(bf, "test_sub.beagle.gz", sep = "\t", col.names = TRUE, row.names = FALSE)
  man3 <- plot_manhattan(list(res, res2$SNP_pvalues), SNP_metadata = "test_sub.beagle.gz")
  expect_equal(dplyr::arrange(man2@data[,.(CHROM, POS, start, cum.bp, nlogp, Env)], CHROM, POS, Env),
               dplyr::arrange(man3@data[,.(CHROM, POS, start, cum.bp, nlogp, Env)], CHROM, POS, Env)) # same data means same plot.


  file.remove("test_sub.beagle.gz")
  unlink("files", recursive = TRUE)
  unlink("plot", recursive = TRUE)
  unlink("stan", recursive = TRUE)
})




