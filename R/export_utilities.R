#' Generate an example SLURM array job script and handlers to run
#' \code{\link{mechGEA}} on a SLURM HPC.
#'
#' Please note: the resulting shell script, located in the \code{output_dir}
#' named \code{example_slurm_array_job.sh} should be edited to fit the user's
#' SLURM configuration. See notes at the top of the file!
#'
#' @param beagle_files vector of file paths to \code{.beagle}/\code{.beagle.gz}
#'   files to run. These will be formatted to phred files for running by the
#'   generated SLURM script for more effective parallelization.
#' @param sample_metadata either a named list of data.frames containing the
#'   needed metadata in the correct format expected by \code{\link{mechGEA}}
#'   (see \code{\link{example_gl}}, although the \code{$likelihoods} object is
#'   not needed.) or an input \code{\link[sf:sf]{sf}} object like
#'   \code{\link{sample_metadata}}.
#' @param output_dir character, default ".". A path to the director where final
#'   \code{\link{mechGEA}} outputs will be stored in named sub-directories for
#'   each file in \code{beagle_files}.
#' @param conda_env character, default mechGEA. A named \code{conda} environment
#'   to load at the start of the script. Note that if a conda module etc. needs
#'   to be loaded for this env to load, the resulting SLURM array-job shell
#'   script must be modified. See \code{\link{setup_mechGEA_env}} for help
#'   setting up a working \code{conda} environment.
#' @param control_file_prefix character, default "./control/". A prefix indicating
#'   where the control files (including meta-data and the shell script) will be
#'   located. Either a directory ending in "/" or a file basename prefix.
#' @param env_mean_method character, default "unweighted". Weighting scheme for
#'   mean envs if a \code{\link{sf}} object is provided to \code{sample_metadata}.
#'   See \code{\link{prepare_metadata}} for details and options.
#' @param prune_r2 numeric or FALSE, default 0.7. If numeric, the r2 by which
#'   to prune environmental variables to ensure all are uncorrelated. Scaling
#'   is always set to TRUE; to scale without pruning set prune_r2 to 1.
#' @param idp integer, default 1. Inverse-distance weighting provided if using
#'   weighted env means. See \code{\link{prepare_metadata}} for details and
#'   options.
#' @param migration_results file path or NULL, default NULL. If provided,
#'   \code{\link{mechGEA}} will use migration results in the list format
#'   produced by \code{\link{estimate_migration}}, saved as an RDS with
#'   \code{\link{saveRDS}} rather than re-running with the provided data. This
#'   is useful if you wish to assess migration once with a set of loci and then
#'   use that assessment for other genomic chunks (see
#'   \code{stop_after_migration} in \code{\link{mechGEA}} to estimate migration
#'   on one chunk and then stop the model).
#' @param memory character, default "100G". Amount of memory to allocate to
#'   each job. See SLURM "--mem" option details for format.
#' @param job_time character, default "1-24:00:00" (two days). Amount of time
#'   to allocate to each job. See SLURM "-t" option details for format.
#' @param cores integer, default 40. Number of cores to allocate for each job.
#'   See SLURM "-n" option details.
#' @param ... additional arguments passed to \code{\link{mechGEA}}.
#'
#' @export
#' @author William Hemstrom
#'
#' @examples
#' # first, chunk a beagle file
#' beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
#' chunks <- get_chunks(beagle, loci_per_chunk = 150)
#'
#' # next, setup a conda env
#' \dontrun{
#' setup_mechGEA_env("mechGEA")
#'
#' # then, create the slurm job
#' slurm_script <- example_slurm_job(
#'   beagle_files = chunks,
#'   sample_metadata = sample_metadata, # the example metadata sf object
#'   conda_env = "mechGEA",
#'   control_file_prefix = "./control/",
#'   env_mean_method = "unweighted", # unweighted env averages
#'   prune_r2 = 0.7 # prune envs
#'   )
#'
#' # lastly, edit the resulting slurm script with cluster/resource specific
#' # needs and submit the job.
#'}
example_slurm_job <- function(beagle_files,
                              sample_metadata,
                              output_dir = ".",
                              conda_env = "mechGEA",
                              control_file_prefix = "./control/",
                              env_mean_method = "unweighted",
                              prune_r2 = .7,
                              idp = 1,
                              migration_results = NULL,
                              memory = "100G",
                              job_time = "1-24:00:00",
                              cores = "40",
                              ...){

  #NOTE: set the resources from arg too?


  beagle_files <- normalizePath(beagle_files)
  bad_beagles <- !file.exists(beagle_files)
  if(any(bad_beagles)){
    stop(paste0("Some beagle files not found:\n", paste0(bad_beagles, collapse = "\n\n")))
  }

  nfiles <- length(beagle_files)

  # create and sanitize output dir
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }

  # sanitize control file prefix (expects trailing / IF A DIRECTORY)
  add_trail <-  grepl("[/\\\\]$", control_file_prefix)
  control_file_prefix <- normalizePath(control_file_prefix, mustWork = FALSE)
  if(add_trail) control_file_prefix <- paste0(control_file_prefix, .Platform$file.sep)

  ## make base dir if needed
  bd <- dirname(control_file_prefix)
  if(!dir.exists(bd)) dir.create(bd)

  # if a dir, make it if needed
  if(grepl("[/\\\\]$", control_file_prefix)){
    if(!dir.exists(control_file_prefix)){
      dir.create(control_file_prefix)
    }
  }

  # save a beagle list
  bfp <- normalizePath(paste0(control_file_prefix, "beaglelist.txt"), mustWork = FALSE)
  writeLines(beagle_files, bfp, sep = "\n")
  beagle_files <- bfp

  # format and save the metadata
  if(is.list(sample_metadata) & !methods::is(sample_metadata, "sf")){
    if(all(c("patch_info", "patches_xy", "env") %in% names(sample_metadata))){

      if(!isFALSE(prune_r2)){
        sample_metadata$envs <- prune_envs(sample_metadata$env, r2 = prune_r2)
      }

      smpath <- normalizePath(paste0(control_file_prefix, "sample_metadata.RDS"), mustWork = FALSE)
      saveRDS(sample_metadata[c("patch_info", "patches_xy", "env")], smpath)
      sample_metadata <- smpath

    } else{
      stop(paste0("'sample_metadata' is a list, but does not contain needed parts:",
                  paste0(c("patch_info", "patches_xy", "env"), collapse = ", "), "\n"))
    }
  } else {
    sample_metadata <- prepare_metadata(sample_metadata, env_mean_method, idp = idp)
    if(!isFALSE(prune_r2)){
      sample_metadata$env <- prune_envs(sample_metadata$env, r2 = prune_r2)
    }
    smpath <- normalizePath(paste0(control_file_prefix, "sample_metadata.RDS"), mustWork = FALSE)
    saveRDS(sample_metadata, smpath)
    sample_metadata <- smpath
  }

  # check the migration data if provided
  if(!is.null(migration_results)){
    migration_results <- normalizePath(migration_results, mustWork = FALSE) # we check it works next line anyway
    if(!file.exists(migration_results)){
      stop(paste0("Could not find pre-existing migration results: ", migration_results, "\n"))
    }
  }

  # read in the shell script
  path <- system.file("extdata", "example_slurm_array_job.sh", package = "mechGEA")
  sh <- paste0(control_file_prefix, "example_slurm_array_job.sh")
  sh <- normalizePath(sh, mustWork = FALSE)
  file.copy(path, sh, overwrite = TRUE)


  # generate the parameter file
  args <- list(...)
  control_parmfile <- paste0(control_file_prefix, "parmfile.RDS")
  normalizePath(control_parmfile, mustWork = FALSE)
  saveRDS(args, control_parmfile)

  # copy over the Rscript
  path <- system.file("extdata", "example_array_rscript.R", package = "mechGEA")
  control_rscript <- paste0(control_file_prefix, "example_array_rscript.R")
  control_rscript <- normalizePath(control_rscript, mustWork = FALSE)
  file.copy(path, control_rscript, overwrite = TRUE)


  # substitute in the provided data
  xfun::gsub_file(sh, "BEAG_REPLACE", beagle_files)
  xfun::gsub_file(sh, "META_REPLACE", sample_metadata)
  xfun::gsub_file(sh, "PARAM_REPLACE", control_parmfile)
  xfun::gsub_file(sh, "DIR_REPLACE", output_dir)
  xfun::gsub_file(sh, "ENV_REPLACE", conda_env)
  xfun::gsub_file(sh, "RSCRIPT_REPLACE", control_rscript)
  xfun::gsub_file(sh, "ARRAY_REPLACE", as.character(nfiles))
  xfun::gsub_file(sh, "MIG_RES_REPLACE", ifelse(!is.null(migration_results), migration_results, "FALSE"))
  xfun::gsub_file(sh, "TIME_REPLACE", job_time)
  xfun::gsub_file(sh, "CORES_REPLACE", cores)
  xfun::gsub_file(sh, "MEM_REPLACE", memory)

  return(sh)
}


#' Prune correlated (environmental) variables.
#'
#' Conducts pruning of envirnonmental values according to correlation by
#' removing the variable correlated with the most other variables iteratively
#' until no correlated pairs remain. Ties removed at random.
#'
#' @param envs numeric matrix/data.frame/data.table of
#' @param r2 numeric, default 0.7. Threshold to use for prunning; correlations
#'   (\code{\link{cor}}) values higher than \code{r2} will result in the removal
#'   of envs.
#' @param scale logical, default TRUE. If TRUE, environmental vairables will
#'   be scaled to mean 0, sd 1. This is recommended for subsequent model fitting
#'   by \code{\link{mechGEA}}.
#' @param seed integer, default a random number between 1 and 10,000. Seed to
#'   use for breaking ties.
#'
#' @export
#' @author William Hemstrom
prune_envs <- function(envs, r2 = .7, scale = TRUE, seed = floor(stats::runif(1, 1, 10000))){
  #============sanity=========
  dreamerr::check_set_value(r2, "scalar NONA numeric GE{0} LT{1}")
  dreamerr::check_set_value(scale, "scalar NONA logical")
  dreamerr::check_set_value(seed, "scalar NONA integer")
  dreamerr::check_set_value(envs, "data.frame | matrix")


  ABcor <- A <- B <- ..keeps <- NULL

  #===========run=============
  set.seed(seed)
  envs <- as.data.table(envs)

  pcol <-FALSE
  if(colnames(envs)[1] == "patch"){
    pcol <- TRUE
    PID <- envs$patch
    envs <- envs[,-1]
  }

  if(scale){
    envs[, (names(envs)) := lapply(.SD, function(x) as.numeric(scale(x)))]
  }

  # purge highly correlated envs
  tcor <- stats::cor(as.matrix(envs))
  tcor <- abs(tcor)
  tcor[lower.tri(tcor, diag = TRUE)] <- 0
  bad_pairs <- tcor >= r2
  bad_pairs <- data.table(ABcor = tcor[bad_pairs],
                          A = colnames(tcor)[col(tcor)[bad_pairs]],
                          B = rownames(tcor)[row(tcor)[bad_pairs]])
  bad_pairs <- bad_pairs[order(ABcor, decreasing = TRUE),]
  all_purged <- character()
  while(nrow(bad_pairs) > 0){
    worst <- table(c(bad_pairs$A, bad_pairs$B))
    worst <- worst[worst == max(worst)]
    purged <- stats::rmultinom(1, 1, prob = rep(1/length(worst), length(worst)))
    purged <- names(worst)[which(purged == 1)]

    bad_pairs <- bad_pairs[A != purged & B != purged,]
    all_purged <- c(all_purged, purged)
  }

  keeps <- which(!colnames(envs) %in% all_purged)
  envs <- .fix..call(envs[,..keeps])

  if(pcol){
    envs <- cbind(patch = PID, as.data.frame(envs))
  } else {
    envs <- as.data.frame(envs)
  }

  return(envs)
}
