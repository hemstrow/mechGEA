test_that("gea plots",{

  skip_if_not_installed("cmdstanr")
  skip_on_cran()

  beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
  phred <- beagle2pl(beagle, lines = 1:150)

  input_data <- prepare_metadata(sample_metadata)

  res <- mechGEA(phred,
                 patch_info = input_data$patch_info,
                 patches_xy = input_data$patches_xy,
                 migration_loci = 50,
                 env = input_data$env,
                 iter_sampling_gea = 100,
                 iter_warmup_gea = 200,
                 iter_sampling_mig = 100,
                 iter_warmup_mig = 200)

  p <- plot_gea(stan_data = res$stan_and_migration_data$stan_data,
                migration_results = res$migration_results,
                beta_draw = res$GEA_results$beta_draw,
                beta_summary = res$GEA_results$beta_summary,
                rho_kappa_mu =  res$GEA_results$rho_kappa_mu,
                env_names = colnames(input_data$env)[-1],
                CI = 0.95)

})
