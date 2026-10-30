test_that("slurm_example", {

  # first, chunk a beagle file
  beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
  chunks <- get_chunks(beagle, loci_per_chunk = 150)

  # then, create the slurm job
  slurm_script <- example_slurm_job(
    beagle_files = chunks,
    sample_metadata = sample_metadata, # the example metadata sf object
    conda_env = "mechGEA",
    control_file_prefix = "./control/",
    env_mean_method = "weighted", # unweighted env averages
    idp = 2,
    output_dir = "./outputs/",
    prune_r2 = 0.7, # prune envs
    iter_sampling_mig = 100,
    iter_warmup_mig = 50,
    max_treedepth_mig = 12
    )

  # see if everything looks OK in the resulting file
  sh <- readLines("control/example_slurm_array_job.sh")

  # correct inputs
  inputs <- sh[grep("Input files", sh):grep("conda activate", sh)]
  inputs <- inputs[-(9:10)]
  ## beagle list
  beag_list <- inputs[grep("beagle", inputs)]
  beag_list <- gsub("beagle_list=", "", beag_list)
  expect_true(file.exists(beag_list))
  flist <- readLines(beag_list)
  expect_true(all(file.exists(flist)))
  ## metadata
  meta <- inputs[grep("metadata", inputs)]
  metaf <- gsub("metadata=", "", meta)
  expect_true(file.exists(metaf))
  meta <- readRDS(metaf)
  expect_equal(names(meta), c("patch_info", "patches_xy", "env"))
  ## parameter file
  par <- inputs[grep("parameter", inputs)]
  parf <- gsub("parameter_file=", "", par)
  expect_true(file.exists(parf))
  par <- readRDS(parf)
  expect_equal(names(par), c("iter_sampling_mig", "iter_warmup_mig", "max_treedepth_mig"))
  ## prefix
  dest_path <- inputs[grep("dir_prefix", inputs)]
  dest_path <- gsub("dir_prefix=", "", dest_path)
  expect_equal(normalizePath(dest_path), normalizePath("./outputs/"))
  ## env
  expect_equal(inputs[grep("env=", inputs)], "env=mechGEA")
  ## rscript
  rfile <- inputs[grep("execution_Rscript", inputs)]
  rfile <- gsub("execution_Rscript=", "", rfile)
  expect_true(file.exists(rfile))
  rlines <- readLines(rfile)
  rlines2 <- readLines(system.file("extdata", "example_array_rscript.R", package = "mechGEA"))
  expect_equal(rlines, rlines2)
  ## no migration provided
  mig <- inputs[grep("mig_res", inputs)]
  expect_equal(mig, "mig_res=FALSE")

  # resources
  ## array
  array <- sh[grep("--array=", sh)]
  array <- as.numeric(gsub("#SBATCH --array=1-", "", array))
  expect_equal(length(flist), array)
  ## time
  time <- sh[grep("#SBATCH -t", sh)]
  time <- gsub("#SBATCH -t ", "", time)
  expect_equal("1-24:00:00", time)
  ## cores
  cores <- sh[grep("#SBATCH -n", sh)]
  cores <- as.numeric(gsub("#SBATCH -n ", "", cores))
  expect_equal(40, cores)
  ## memory
  mem <- sh[grep("#SBATCH --mem=", sh)]
  mem <- gsub("#SBATCH --mem=", "", mem)
  expect_equal("100G", mem)


  # cleanup
  file.remove(rfile)
  file.remove(metaf) # cleanup
  file.remove(parf) # cleanup
  file.remove(flist) # cleanup
  unlink("./control/", recursive = TRUE)



  # once more, but with provided mig as a test
  ## in this case, just provide a BS mig file to check it fills correctly
  BS_mig <- numeric()
  saveRDS(BS_mig, "BS_mig.RDS")


  beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
  chunks <- get_chunks(beagle, loci_per_chunk = 150)

  # then, create the slurm job
  slurm_script <- example_slurm_job(
    beagle_files = chunks,
    sample_metadata = sample_metadata, # the example metadata sf object
    conda_env = "mechGEA",
    control_file_prefix = "./control/",
    env_mean_method = "weighted", # unweighted env averages
    idp = 2,
    output_dir = "./outputs/",
    prune_r2 = 0.7, # prune envs
    iter_sampling_mig = 100,
    iter_warmup_mig = 50,
    max_treedepth_mig = 12,
    migration_results = "BS_mig.RDS"
  )

  sh <- readLines("control/example_slurm_array_job.sh")

  # correct inputs
  inputs <- sh[grep("Input files", sh):grep("conda activate", sh)]
  inputs <- inputs[-(9:10)]
  mig <- inputs[grep("mig_res", inputs)]
  mig <- gsub("mig_res=", "", mig)
  expect_true(file.exists(mig))

  ## once more, but with a mig file that DOESN'T exist
  expect_error(slurm_script <- example_slurm_job(
    beagle_files = chunks,
    sample_metadata = sample_metadata, # the example metadata sf object
    conda_env = "mechGEA",
    control_file_prefix = "./control/",
    env_mean_method = "weighted", # unweighted env averages
    idp = 2,
    output_dir = "./outputs/",
    prune_r2 = 0.7, # prune envs
    iter_sampling_mig = 100,
    iter_warmup_mig = 50,
    max_treedepth_mig = 12,
    migration_results = "BS_mig_NOT_REAL.RDS"
  ), "Could not find pre-existing migration results")


  # cleanup
  file.remove(rfile)
  file.remove(metaf) # cleanup
  file.remove(parf) # cleanup
  file.remove(flist) # cleanup
  unlink("./control/", recursive = TRUE)
  file.remove("BS_mig.RDS")
})


test_that("env setup",{

  skip_on_cran(); skip_on_ci()

  # pulled from snpR internal
  rand_strings <- function (n, length) {
    chrs <- sample(c(LETTERS, letters, 1:9), n * length, replace = TRUE)
    chrs <- split(chrs, rep(1:n, length.out = n * length))
    chrs <- unlist(lapply(chrs, function(x) paste0(x, collapse = "")))
    chrs <- as.character(chrs)
    return(chrs)
  }

  # name for the test env install
  env_name <- rand_strings(1, 5)
  env_exists <- env_name %in% reticulate::conda_list()$name
  while(env_exists){
    env_name <- rand_strings(1, 5)
    env_exists <- env_name %in% reticulate::conda_list()$name
  }

  # test running the env; slow AF.
  setup_mechGEA_env(env_name)

  # run a test in the env
  writeLines(sprintf('
library(mechGEA)

beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
phred  <- beagle2pl(beagle, lines = 1:150)

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

saveRDS(res, "%s")
', "res.RDS"), "test_mechGEA.R")

  reticulate::conda_run2(
    cmd = "Rscript", args = "test_mechGEA.R",
    envname = env_name)

  res <- readRDS("res.RDS")
  expect_true(sum(res$SNP_pvalues$nlogp > 4) == 5)

  # cleanup
  reticulate::conda_remove(envname = env_name)
  file.remove("res.RDS")
  file.remove("test_mechGEA.R")
})
