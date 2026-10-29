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
