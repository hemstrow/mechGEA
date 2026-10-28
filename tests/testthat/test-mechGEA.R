test_that("full pipeline",{

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





  # check the many, many outputs. For now, these are mostly format checks on terminal outputs.
  # NIKUNJ, good place to add here.
  expect_true(sum(res$SNP_pvalues$nlogp > 4) == 5)
  expect_true(nrow(res$GEA_results$beta_summary) == 150*5)
  expect_equal(dim(res$migration_results), c(20, 20))
  expect_true(all(round(rowSums(res$migration_results), 3) == 0))
  expect_equal(dim(res$GEA_results$rho_kappa_mu), c(100*2, 450))
  expect_equal(dim(res$GEA_results$beta_summary), c(150*5, 10))
  expect_equal(dim(res$GEA_results$tau), c(100*2, 1))
  expect_equal(dim(res$GEA_results$mlog_p), c(150, 5))
  expect_equal(dim(res$Neff_results$G_eff_stan), c(150, 20))
  expect_equal(dim(res$Neff_results$G_eff_ana), c(150, 20))
  expect_equal(dim(res$Neff_results$TwoN_eff_stan), c(150, 20))
  expect_equal(dim(res$Neff_results$TwoN_eff_ana), c(150, 20))
  expect_true(cor(as.numeric(res$Neff_results$G_eff_stan), as.numeric(res$Neff_results$G_eff_ana)) > 0.95) # should be pretty correlated
  expect_true(is.numeric(res$Neff_results$Neff_mean))
  expect_true(length(res$Neff_results$Neff_mean) == 1)
  expect_true(res$Neff_results$Neff_mean <= 1)
})
