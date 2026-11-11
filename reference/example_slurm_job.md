# Generate an example SLURM array job script and handlers to run [`mechGEA`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md) on a SLURM HPC.

Please note: the resulting shell script, located in the `output_dir`
named `example_slurm_array_job.sh` should be edited to fit the user's
SLURM configuration. See notes at the top of the file!

## Usage

``` r
example_slurm_job(
  beagle_files,
  sample_metadata,
  output_dir = ".",
  conda_env = "mechGEA",
  control_file_prefix = "./control/",
  env_mean_method = "unweighted",
  prune_r2 = 0.7,
  idp = 1,
  migration_results = NULL,
  memory = "100G",
  job_time = "1-24:00:00",
  cores = "40",
  setup_env = FALSE,
  setup_env_conda_path = "DEFAULT",
  ...
)
```

## Arguments

- beagle_files:

  vector of file paths to `.beagle`/`.beagle.gz` files to run. These
  will be formatted to phred files for running by the generated SLURM
  script for more effective parallelization.

- sample_metadata:

  either a named list of data.frames containing the needed metadata in
  the correct format expected by
  [`mechGEA`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md)
  (see
  [`example_gl`](https://hemstrow.github.io/mechGEA/reference/example_gl.md),
  although the `$likelihoods` object is not needed.) or an input 'sf'
  object like
  [`sample_metadata`](https://hemstrow.github.io/mechGEA/reference/sample_metadata.md).

- output_dir:

  character, default ".". A path to the director where final
  [`mechGEA`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md)
  outputs will be stored in named sub-directories for each file in
  `beagle_files`.

- conda_env:

  character, default mechGEA. A named `conda` environment to load at the
  start of the script. Note that if a conda module etc. needs to be
  loaded for this env to load, the resulting SLURM array-job shell
  script must be modified. See
  [`setup_mechGEA_env`](https://hemstrow.github.io/mechGEA/reference/setup_mechGEA_env.md)
  for help setting up a working `conda` environment. If `conda_env` is
  FALSE, no conda environment will be used.

- control_file_prefix:

  character, default "./control/". A prefix indicating where the control
  files (including meta-data and the shell script) will be located.
  Either a directory ending in "/" or a file basename prefix.

- env_mean_method:

  character, default "unweighted". Weighting scheme for mean envs if a
  'sf' object is provided to `sample_metadata`. See
  [`prepare_metadata`](https://hemstrow.github.io/mechGEA/reference/prepare_metadata.md)
  for details and options.

- prune_r2:

  numeric or FALSE, default 0.7. If numeric, the r2 by which to prune
  environmental variables to ensure all are uncorrelated. Scaling is
  always set to TRUE; to scale without pruning set prune_r2 to 1.

- idp:

  integer, default 1. Inverse-distance weighting provided if using
  weighted env means. See
  [`prepare_metadata`](https://hemstrow.github.io/mechGEA/reference/prepare_metadata.md)
  for details and options.

- migration_results:

  file path or NULL, default NULL. If provided,
  [`mechGEA`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md)
  will use migration results in the list format produced by
  [`estimate_migration`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md),
  saved as an RDS with [`saveRDS`](https://rdrr.io/r/base/readRDS.html)
  rather than re-running with the provided data. This is useful if you
  wish to assess migration once with a set of loci and then use that
  assessment for other genomic chunks (see `stop_after_migration` in
  [`mechGEA`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md)
  to estimate migration on one chunk and then stop the model).

- memory:

  character, default "100G". Amount of memory to allocate to each job.
  See SLURM "–mem" option details for format.

- job_time:

  character, default "1-24:00:00" (two days). Amount of time to allocate
  to each job. See SLURM "-t" option details for format.

- cores:

  integer, default 40. Number of cores to allocate for each job. See
  SLURM "-n" option details.

- setup_env:

  logical, default FALSE. If TRUE, will automatically attempt to setup a
  conda environment with the name given in `conda_env` using
  [`setup_mechGEA_env`](https://hemstrow.github.io/mechGEA/reference/setup_mechGEA_env.md).
  If a conda environment with that name already exists, the user will be
  prompted to overwrite or the function will fail.

- setup_env_conda_path:

  character, default "DEFAULT". Sets the path to the conda executable to
  use; see
  [`setup_mechGEA_env`](https://hemstrow.github.io/mechGEA/reference/setup_mechGEA_env.md).

- ...:

  additional arguments passed to
  [`mechGEA`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md).

## Author

William Hemstrom

## Examples

``` r
# first, chunk a beagle file
beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
chunks <- get_chunks(beagle, loci_per_chunk = 150)
#> Splitting into 10 files with approx. 150 snps per file.

# next, setup a conda env
if (FALSE) { # \dontrun{
setup_mechGEA_env("mechGEA")

# then, create the slurm job
slurm_script <- example_slurm_job(
  beagle_files = chunks,
  sample_metadata = sample_metadata, # the example metadata sf object
  conda_env = "mechGEA",
  control_file_prefix = "./control/",
  env_mean_method = "unweighted", # unweighted env averages
  prune_r2 = 0.7 # prune envs
  )

# lastly, edit the resulting slurm script with cluster/resource specific
# needs and submit the job.
} # }
```
