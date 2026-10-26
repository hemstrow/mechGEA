test_that("multiplication works", {
  expect_error(stop("hi, this test failed"), "this test failed")
})


test_that("full pipeline",{
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

})
