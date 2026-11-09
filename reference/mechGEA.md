# Mechanistic genome-environment association tests.

Core functions for the mechanistic GEA-with-migration pipeline,
producing output files and plots in the working directory (and returned
as an object). For large data (e.g. many SNPs), consider running this on
small `Phred` batches and then combining.

## Usage

``` r
mechGEA(
  Phred,
  patch_info,
  patches_xy,
  env,
  migration_results = NULL,
  N_sim_mig = 100,
  migration_loci = 1000,
  k_sd = 4,
  p0 = 20,
  seed = 42,
  migration_model = "migration_fast",
  max_attempts_mig = 10,
  target_div_pct_mig = 1,
  target_rhat_pct_mig = 1,
  chains_mig = 4,
  iter_sampling_mig = 1000,
  iter_warmup_mig = 1000,
  max_treedepth_mig = 12,
  refresh_mig = 200,
  mig_plot_threshold = 0.4,
  mig_plot_topn = 2,
  mig_plot_binsize = 1,
  chains_gea = 2,
  refresh_gea = 200,
  iter_warmup_gea = 250,
  iter_sampling_gea = 500,
  CI = 0.95,
  env_colors = NULL,
  stop_after_migration = FALSE,
  par = parallel::detectCores() - 2
)

prepare_mig_data(
  Phred,
  patch_info,
  patches_xy,
  env,
  migration_loci,
  N_sim_mig,
  seed = 42
)

estimate_migration(
  migration_data,
  migration_model = "migration_fast",
  max_attempts_mig = 10,
  target_div_pct_mig = 1,
  target_rhat_pct_mig = 1,
  chains_mig = 4,
  iter_sampling_mig = 1000,
  iter_warmup_mig = 1000,
  max_treedepth_mig = 12,
  refresh_mig = 200,
  par = parallel::detectCores() - 2
)

run_stan_GEA_model(
  stan_data,
  migration_results,
  k_sd = 4,
  p0 = 20,
  refresh_gea = 100,
  chains_gea = 2,
  iter_warmup_gea = 250,
  iter_sampling_gea = 500,
  seed = 42,
  par = parallel::detectCores() - 2
)

estimate_neff(stan_data)

beta_summary(beta_draw, nloci, nenv)

plot_migration(
  migration_results,
  migration_data,
  patches_xy,
  mig_plot_threshold = 0.4,
  mig_plot_topn = 2,
  mig_plot_binsize = 1
)

plot_gea(
  stan_data,
  migration_results,
  beta_draw,
  beta_summary,
  rho_kappa_mu,
  env_names,
  CI = 0.95,
  env_colors = NULL
)
```

## Arguments

- Phred:

  Phred-scaled genotype liklelihoods, as generated from a `.beagle` or
  `.beagle.gz` file by
  [`beagle2pl`](https://hemstrow.github.io/mechGEA/reference/beagle2pl.md).
  See
  [`example_gl`](https://hemstrow.github.io/mechGEA/reference/example_gl.md)`$likelihoods`
  for an example.

- patch_info:

  [`data.frame`](https://rdrr.io/r/base/data.frame.html) or similar.
  Patch information for each sample in a two-column data.frame with
  column names `ind_id` and `patch`. Important note: samples *MUST* be
  sorted in blocks by patch (all samples in one patch must be in
  blocks). Data fitting this format can be produced by
  [`prepare_metadata`](https://hemstrow.github.io/mechGEA/reference/prepare_metadata.md).

- patches_xy:

  [`data.frame`](https://rdrr.io/r/base/data.frame.html) or similar.
  Spatial coordinates for each patch, with three columns named `patch`,
  `lat`, and `lon`. Must be ordered identically to order of patch
  appearance in `patch_info`. Data fitting this format can be produced
  by
  [`prepare_metadata`](https://hemstrow.github.io/mechGEA/reference/prepare_metadata.md).

- env:

  [`data.frame`](https://rdrr.io/r/base/data.frame.html) or similar.
  Uncorrelated, scaled environmental variables for each patch, with N +
  1 columns. The first column must be named `patch` and contain patch
  names. The additional N columns can have any name and can contain
  Z-scaled env data (mean = 0, sd = 1). Patches must be ordered
  identically to order of patch appearance in `patch_info`. Data fitting
  this format can be produced by
  [`prepare_metadata`](https://hemstrow.github.io/mechGEA/reference/prepare_metadata.md)
  and uncorrelated envs can be selected using
  [`prune_envs`](https://hemstrow.github.io/mechGEA/reference/prune_envs.md).

- migration_results:

  list or NULL, default NULL. If provided, `mechGEA` will use migration
  results in the list format produced by `estimate_migration` rather
  than re-running with the provided data. This is useful if you wish to
  assess migration once with a set of loci and then use that assessment
  for other genomic chunks (see `stop_after_migration` to estimate
  migration and then stop the model).

- N_sim_mig:

  integer, default 100. Number of times genotypes are re-sampled when
  estimating the diversity for migration. Ignored if `migration_results`
  are provided.

- migration_loci:

  integer, default 1000. Number of loci to use from `Phred` to use for
  estimating migration. Must be equal to or less than the total number
  of provided loci. Ignored if `migration_results` are provided; if less
  than 1000, will be set to 1000 or the number of loci in the dataset.

- k_sd:

  numeric, default 4. NIKUNJ HERE PLEASE!

- p0:

  numeric, default 20. NIKUNJ HERE PLEASE!

- seed:

  integer, default 42. Seed used for all model steps.

- migration_model:

  character, default "migration_fast". Choice of model for estimating
  migration. Options are "migration_fast" or "migration_M3b". See
  details. NIKUNJ HERE PLEASE!

- max_attempts_mig:

  integer, default 10. Number of times to try fitting the migration
  model. If the model does not converge after `max_attempts_mig` tries,
  the best fit model will be used.

- target_div_pct_mig:

  numeric, default 1. Maximum percentage of divergences allowed during
  migration model fitting.

- target_rhat_pct_mig:

  numeric, default 1. Maximum percentage of rhat \> 1.05 allowed during
  fitting.

- chains_mig:

  integer, default 4. Number of chains to run for the Bayesian migration
  model fitting.

- iter_sampling_mig:

  integer, default 1000. Number of sampling iterations to run the
  Bayesian migration model fitting.

- iter_warmup_mig:

  integer, default 1000. Number of warm-up/burn-in iterations to run the
  Bayesian migration model fitting.

- max_treedepth_mig:

  integer, default 12. Maximum tree-depth to use during Bayesian
  migration model fitting.

- refresh_mig:

  integer, default 200. Controls how often model progress is reported
  during Bayesian migration model fitting.

- mig_plot_threshold:

  numeric, default 0.4. For plotting migration model results, controls
  the migration strength threshold above which migration connections
  will be included in networks.

- mig_plot_topn:

  integer, default 2. Number of migratory connections to consider for
  each pixel when plotting migration networks. NIKUNJ HERE PLEASE!

- mig_plot_binsize:

  integer, default 1. NIKUNJ HERE PLEASE!

- chains_gea:

  integer, default 2. Number of chains to run for the GEA model.

- refresh_gea:

  integer, default 200. Controls how often model progress is reported
  during GEA model fitting.

- iter_warmup_gea:

  integer, default 200. Number of warm-up/burn-in iterations to run the
  during GEA model fitting.

- iter_sampling_gea:

  integer, default 250. Number of sampling iterations to run the during
  GEA model fitting.

- CI:

  numeric, default 0.95. Confidence interval to plot for the GEA model.

- env_colors:

  character or NULL, default NULL. Vector of colors for environmental
  variables.

- stop_after_migration:

  logical, default FALSE. If TRUE, stop the model after running
  `estimate_migration`. This is useful if you wish to assess migration
  once with a set of loci and then use that assessment for other genomic
  chunks (see `migration_results` to pick up the model with the results
  in other chunks.).

- par:

  integer, default `parallel::detectCores() - 2`, two less than the
  available cores. Maximum number of cores to allocate for model
  fitting; the number of cores per chain will be `par`/`N_chains`.

- migration_data:

  list containing data for migration model, in the format produced by
  `prepare_mig_data``$migration_data`.

- stan_data:

  list containing data for the stan GEA model fitting, in the format
  produced by `prepare_mig_data``$stan_data`.

- beta_draw:

  beta parameter sampling/draw data from a fit GEA model, as produced by
  `run_stan_GEA_model``$beta_draw`.

- nloci:

  integer, number of loci in the source data for the stan model.

- nenv:

  integer, number of environmental variables in the source data for the
  stan model.

- beta_summary:

  beta parameter summary data from a fit GEA model, as produced by
  `beta_summary``$beta_summary`.

- rho_kappa_mu:

  parameter fit data from a GEA model, as produced by
  `run_stan_GEA_model``$rho_kapp_mu`.

- env_names:

  character, vector of environmental variable names.

## Value

A nested list containing model outputs. Different modules may return
different portions of the following list; `mechGEA` will return every
part. NIKUNJ, DETAIL WOULD BE GREAT!

- SNP_pvalues: A data.table containing associations (`nlogp`,
  \\-\log\_{e}\\ transformed p-values) for each SNP (numbered 1:n by
  input order `Phred`) compared to each environmental variable.

- GEA_results: A list containing summaries of the GEA model fits.

- migration_results: A n by n matrix where n is the number of patches;
  estimated migration rates between each patch.

- Neff_results: A list containing the effective sample size estimates
  for each locus vs each environmental variable.

- migration_plot: Plots of isolation-by-distance and migration rates
  between patches.

- GEA_plots: Manhattan plot, model fit visualization, and
  migration-selection balances for each significant correlation from the
  GEA.

- stan_and_migration_data: Formatted data which was used by the GEA and
  migration models.

## Functions

- `mechGEA()`: The full mechanistic GEA model pipeline.

- `prepare_mig_data()`: Prepare migration data (and stan) data for
  modeling. Consider running with `mechGEA` for the full pipeline.

- `estimate_migration()`: Run the migration estimation model. Consider
  running with `mechGEA` for the full pipeline.

- `run_stan_GEA_model()`: Run the Bayesian mechanistic GEA model.
  Consider running with `mechGEA` for the full pipeline.

- `estimate_neff()`: Estimate effective sample sizes for data for
  `run_stan_GEA_model`. Consider running with `mechGEA` for the full
  pipeline.

- `beta_summary()`: Calculate parameter summaries from a
  `run_stan_GEA_model`. Consider running with `mechGEA` for the full
  pipeline.

- `plot_migration()`: Plot migration model results. Consider running
  with `mechGEA` for the full pipeline.

- `plot_gea()`: Plot GEA model results. Consider running with `mechGEA`
  for the full pipeline.

## References

Goel, N., et al. (2026). Identifying adaptive variation in spatially
structured populations using low-coverage whole-genome sequencing data.
*Molecular Biology and Evolution*, 43(10), msag255.
https://doi.org/10.1093/molbev/msag255

## Author

William Hemstrom, Nikunj Goel

## Examples

``` r
# a basic run of the full pipeline
## example beagle file:
beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
## convert the formats of the data and metadata
phred <- beagle2pl(beagle, lines = 1:150)
input_data <- prepare_metadata(sample_metadata)
## run the pipeline
if (FALSE) { # \dontrun{
res <- mechGEA(phred,
               patch_info = input_data$patch_info,
               patches_xy = input_data$patches_xy,
               migration_loci = 50,
               env = input_data$env,
               iter_sampling_gea = 100,
               iter_warmup_gea = 200,
               iter_sampling_mig = 100,
               iter_warmup_mig = 200)
} # }

# Run only up until the migration calculation and plot results
if (FALSE) { # \dontrun{
res_mig <- mechGEA(phred,
                   patch_info = input_data$patch_info,
                   patches_xy = input_data$patches_xy,
                   migration_loci = 50,
                   env = input_data$env,
                   iter_sampling_mig = 100,
                   iter_warmup_mig = 200,
                   stop_after_migration = TRUE)
p <- plot_migration(migration_results = res_mig$migration_results,
       migration_data = res_mig$stan_and_migration_data$migration_data,
       patches_xy = input_data$patches_xy,
       mig_plot_threshold = 0.4,
       mig_plot_topn = 2,
       mig_plot_binsize = 1)
} # }

# using a pre-computed migration result (great if running a migration once
# and applying across many chunks)
if (FALSE) { # \dontrun{
res_ <- mechGEA(phred,
               patch_info = input_data$patch_info,
               patches_xy = input_data$patches_xy,
               migration_loci = 50,
               env = input_data$env,
               iter_sampling_gea = 100,
               iter_warmup_gea = 200,
               migration_results = res_mig$migration_results)
} # }

```
