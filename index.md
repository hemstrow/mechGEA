# mechGEA

# Background

**mechGEA** is a mechanistic genotype–environment association (GEA)
pipeline for identifying loci involved in local adaptation, designed to
work directly with low-coverage sequencing data.

## Why a mechanistic GEA?

Most GEA methods have two limitations:

1.  **They ignore uncertainty in genotyping and cannot have missing
    data.** Allele frequencies are estimated from a finite, often
    incomplete sample of genotypes, and with low-coverage data the
    genotypes themselves are uncertain. Treating these estimates as
    known can bias downstream inference. Similarly, they typically also
    do not allow missing data, requiring imputation or a similar
    approach prior to use.
2.  **They are phenomenological.** Regression-based approaches
    (e.g. GLMs) assume a particular shape for the relationship between
    environment and allele frequency, but are only loosely connected to
    evolutionary theory.

mechGEA instead builds the analysis on an explicit model of evolution
and incorporates uncertainty into the model directly via data model.

# Installation

mechGEA fits models with Stan via `cmdstanr` and `rstan`. Installing
these dependencies manually can be difficult, in part because they must
be installed in a particular order under different C++ compiler
versions. The usual installation tools (`install.packages` or the
`remotes` and `pak` packages) cannot handle this on their own;
installing with these tools alone will only produce a partial-install.

We therefore **strongly recommend installing mechGEA into a dedicated
conda environment** using a two-step process via the function
[`setup_mechGEA_env()`](https://hemstrow.github.io/mechGEA/reference/setup_mechGEA_env.md)
included in mechGEA.

> mechGEA is currently tested and maintained for **Linux** only. Support
> for Windows and macOS is limited.

**1. Initial install of mechGEA** from GitHub, using `remotes` or `pak`:

``` r

install.packages("remotes") # skip if already installed!
remotes::install_github("hemstrow/mechGEA")
```

This will install the package, but the core functions will be turned off
unless mechGEA detects a complete install.

**2. Build the complete environment** with
[`setup_mechGEA_env()`](https://hemstrow.github.io/mechGEA/reference/setup_mechGEA_env.md):

``` r

mechGEA::setup_mechGEA_env() # creates a conda environment called "mechGEA"
```

This builds a self-contained conda environment (using the `reticulate`
package) with everything mechGEA needs: R, C/C++/Fortran compilers,
CmdStan, `cmdstanr`, `rstan`, `sf`, `terra`, all other dependencies, and
mechGEA itself. No system libraries or administrator rights are
required. - If conda is not found, Miniconda is installed automatically.
To use an existing conda installation (other than the default), provide
a path to the conda executable:
`setup_mechGEA_env(conda_path = "/path/to/conda")`. - Use `name` to
choose a different environment name. The function stops if an
environment with that name already exists. - Setup compiles several
packages and can take a long time. At the end, check that the printed
CmdStan path points inside the conda environment.

**3. Activate the environment** in bash to run analyses, then either run
an R script or start up an R session:

``` bash
conda activate mechGEA
Rscript path/to/my_analysis.R
```

# Quick-start Guide

Here, we will walk through a basic analysis of the example data included
in this package (`example_gl` and `sample_metadata`): 200 individuals
sampled across 20 patches, genotyped at 1,500 loci, with five
environmental variables. We will look at the input formats, prepare
patch metadata, estimate migration among patches, fit the mechanistic
GEA model, and summarise and plot the results. We first run the full
pipeline with a single call to
[`mechGEA()`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md),
then step through each stage individually.

``` r

library(mechGEA)
#> Loading required package: data.table
#> Warning: package 'data.table' was built under R version 4.5.3
#> mechGEA version: 0.0.0.9000 loaded.
#> ================================================================================
#> This is a development build of mechGEA.
#> While this version is passing all tests, there may still be bugs present.
#> Please report any issues to the github issues page (https://github.com/hemstrow/mechGEA/issues)!
#> Check out the NEWS page on github for update notes (https://github.com/hemstrow/mechGEA/blob/dev/NEWS.md)
#> ================================================================================
```

## Data import

mechGEA needs four inputs:

| Input | Description |
|----|----|
| `Phred` | Phred-scaled genotype likelihoods, one row per locus and one column per individual; genotypes comma separated with missing data marked “.” |
| `patch_info` | Patch identification for each individual |
| `patches_xy` | Coordinates of each patch |
| `env` | Z-scaled (mean 0, sd 1), uncorrelated environmental variables for each patch |

For the example data, all four are stored in the `example_gl` list:

> **Important:** Individuals in `sample_metadata` must be in the **same
> order** as the columns (individuals) of the genotype likelihoods, and
> all individuals from the same patch must be in **consecutive rows**.
> Patches must in `patches_xy` and `env` must also be sorted by
> order-of-appearance in `patch_info`.

``` r

str(example_gl, max.level = 1)
#> List of 4
#>  $ likelihoods: chr [1:1500, 1:200] "46,6,0" "138,18,0" "." "69,9,0" ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ patch_info :'data.frame': 200 obs. of  2 variables:
#>  $ patches_xy :'data.frame': 20 obs. of  3 variables:
#>  $ env        :'data.frame': 20 obs. of  6 variables:
```

### Generating correctly formatted metadata from an `sf` object.

If your sample metadata is stored per individual in an `sf` object with
one row per individual and environmental values attached, as in the
`sample_metadata` example data,
the[`prepare_metadata()`](https://hemstrow.github.io/mechGEA/reference/prepare_metadata.md)
function can be used to convert it into `patch_info`, `patches_xy`, and
`env` objects:

``` r

head(sample_metadata)
#>      patch  ind_id        env1      env2       env3     env4     env5
#> 1 patch_01 ind_001 -0.03212425 0.2021607 -0.6491746 1.558552 1.141542
#> 2 patch_01 ind_002 -0.03212425 0.2021607 -0.6491746 1.558552 1.141542
#> 3 patch_01 ind_003 -0.03212425 0.2021607 -0.6491746 1.558552 1.141542
#> 4 patch_01 ind_004 -0.03212425 0.2021607 -0.6491746 1.558552 1.141542
#> 5 patch_01 ind_005 -0.03212425 0.2021607 -0.6491746 1.558552 1.141542
#> 6 patch_01 ind_006 -0.03212425 0.2021607 -0.6491746 1.558552 1.141542
#>               geometry
#> 1 -105.65012, 46.04031
#> 2 -105.65012, 46.04031
#> 3 -105.65012, 46.04031
#> 4 -105.65012, 46.04031
#> 5 -105.65012, 46.04031
#> 6 -105.65012, 46.04031
```

``` r

meta <- prepare_metadata(sample_metadata)
```

It is important to note that the GEA model assumes environmental
variables are uncorrelated and Z-scaled. The
[`prune_envs()`](https://hemstrow.github.io/mechGEA/reference/prune_envs.md)
function can be used to prune highly correlated environmental variables
according to a correlation threshold by sequentially removing the
variables correlated with the highest number of other variables.
[`prune_envs()`](https://hemstrow.github.io/mechGEA/reference/prune_envs.md)
also Z-scales data.

``` r

# removes one env from any pair with a correlation coefficient above 0.7
# (or below -0.7).
meta$env <- prune_envs(meta$env, r2 = 0.7)
```

The resulting, correctly formatted data should be essentially identical
to that in the `example_gl` data set, albeit without genotypes.

### Patch metadata details:

`patch_info` is a two-column data.frame with columns `ind_id` and
`patch`. Its rows must be in the same order as the columns of `Phred`,
and **individuals must be sorted in blocks by patch**: all individuals
from one patch must be in a block (sorted next to each other).

``` r

head(example_gl$patch_info)
#>    ind_id    patch
#> 1 ind_001 patch_01
#> 2 ind_002 patch_01
#> 3 ind_003 patch_01
#> 4 ind_004 patch_01
#> 5 ind_005 patch_01
#> 6 ind_006 patch_01
```

`patches_xy` gives coordinates for each patch, with columns `patch`,
`lat` and `lon`. `env` has a `patch` column followed by one column per
environmental variable, each Z-scaled (mean = 0, sd = 1). Lat/lon
coordinates are recommended to be in “WGS 84 (EPSG:4326)” coordinate
reference for automatic plotting, but this is not required. Both must
list patches in the same order in which patches first appear in
`patch_info`.

``` r

head(example_gl$patches_xy)
#>      patch      lat       lon
#> 1 patch_01 46.04031 -105.6501
#> 2 patch_02 38.38710 -105.3495
#> 3 patch_03 46.88892 -114.1371
#> 4 patch_04 46.46668 -106.7890
#> 5 patch_05 37.82438 -109.3364
#> 6 patch_06 42.14212 -110.9922
head(example_gl$env)
#>      patch        env1       env2        env3         env4        env5
#> 1 patch_01 -0.03212425  0.2021607 -0.64917455  1.558552048  1.14154181
#> 2 patch_02 -1.36106996 -0.3518810  0.03029178  0.191683584  0.99113613
#> 3 patch_03  0.08928409  0.7416528  0.51805408  0.007063579 -0.98067315
#> 4 patch_04  1.33885570 -0.7091368  1.52398967 -0.220972056  1.76498136
#> 5 patch_05  1.95212686 -1.3359885 -1.09200225 -1.390288050 -0.65674777
#> 6 patch_06 -0.14371808  0.4237745  1.40445240  0.577386601  0.08682305
```

A quick sanity check before running anything:

``` r

patch_order <- unique(example_gl$patch_info$patch)
identical(patch_order, example_gl$patches_xy$patch)
#> [1] TRUE
identical(patch_order, example_gl$env$patch)
#> [1] TRUE
identical(colnames(example_gl$likelihoods), example_gl$patch_info$ind_id)
#> [1] TRUE
```

All three should be `TRUE`. Note that these checks are also run
automatically.

## Genotype likelihoods

Rather than called genotypes, mechGEA uses genotype likelihoods, so
uncertainty from low read depth is carried through the model. Each cell
holds three comma-separated Phred-scaled likelihoods for the three
possible genotypes (homozygous reference, heterozygote, homozygous
alternate). Phred-scaling follows the equation: with the lowest
Phred-scaled genotype (the most likely genotype) assigned a zero and
everything else scaled to that value.

``` r

dim(example_gl$likelihoods)          # loci x individuals
#> [1] 1500  200
example_gl$likelihoods[1:5, 1:5]
#>            ind_001    ind_002    ind_003    ind_004   ind_005   
#> locus_0001 "46,6,0"   "0,9,69"   "115,15,0" "8,0,77"  "0,9,69"  
#> locus_0002 "138,18,0" "14,0,37"  "14,0,37"  "46,6,0"  "46,6,0"  
#> locus_0003 "."        "161,21,0" "161,21,0" "37,0,14" "31,0,54" 
#> locus_0004 "69,9,0"   "138,18,0" "138,18,0" "2,0,117" "46,6,0"  
#> locus_0005 "54,0,31"  "92,12,0"  "0,18,138" "92,12,0" "184,24,0"
```

Genotype likelihoods in BEAGLE format (for example, `.beagle.gz` output
from ANGSD) can be converted to this format with
[`beagle2pl()`](https://hemstrow.github.io/mechGEA/reference/beagle2pl.md).
Each cell holds three Phred-scaled genotype likelihoods (AA, AB, BB),
(again, rescaled so the most likely genotype is 0). Correspondingly, a
phred-likelihood of 10 means that genotype is 10× less likely than the
likelihood with a score of 0. “.” marks missing data.

An example BEAGLE file is included with the package:

``` r

bfile <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
Phred <- beagle2pl(bfile)
dim(Phred)
#> [1] 1500  200
```

[`beagle2pl()`](https://hemstrow.github.io/mechGEA/reference/beagle2pl.md)
does not keep sample names from the beagle file by default (since many
programs, like ANGSD, name individuals “Ind0” through “Ind(n-1)”). The
[`beagle2pl()`](https://hemstrow.github.io/mechGEA/reference/beagle2pl.md)
function allows for sample names to be provided via the `ind_names`
argument.

**Note: this assumes `sample_metadata` is in the same order as the BAM
list used to create the BEAGLE file.** If it is not, provide `ind_names`
in the BAM list order and use the `ind_order` argument to resort
individuals on read-in. This can also be used to re-sort individuals
which are not in blocks by patch.

``` r

colnames(Phred) <- sample_metadata$ind_id
Phred[1:5, 1:5]
#>      ind_001    ind_002    ind_003    ind_004   ind_005   
#> [1,] "46,6,0"   "0,9,69"   "115,15,0" "8,0,77"  "0,9,69"  
#> [2,] "138,18,0" "14,0,37"  "14,0,37"  "46,6,0"  "46,6,0"  
#> [3,] "."        "161,21,0" "161,21,0" "37,0,14" "31,0,54" 
#> [4,] "69,9,0"   "138,18,0" "138,18,0" "2,0,117" "46,6,0"  
#> [5,] "54,0,31"  "92,12,0"  "0,18,138" "92,12,0" "184,24,0"
```

# Running the full pipeline

[`mechGEA()`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md)
runs every step of the analysis in one call:

1.  **Data preparation**: validates the inputs, converts Phred scores to
    log-likelihoods, and builds the data for both models.
2.  **Migration estimation**: fits the Bayesian migration model on a
    random subset of loci.
3.  **GEA model fitting**: fits the mechanistic GEA model to all loci,
    accounting for the estimated migration.
4.  **Summaries and plots**: identifies adaptive (locus, environment)
    pairs and plots the results.
5.  **Effective sample size**: estimates effective sample sizes for each
    locus/ environmental variable pair.

Results are returned as a list and also written to `files/` and `plot/`
in the **current working directory**, so run it in a dedicated folder.

The settings below are deliberately small so that the example finishes
quickly (about an hour on a 22-core server). They are **not** suitable
for real analyses; use the defaults, or larger values, for those.

``` r


res <- mechGEA(Phred      = example_gl$likelihoods,
               patch_info = example_gl$patch_info,
               patches_xy = example_gl$patches_xy,
               env        = example_gl$env,
               migration_loci    = 1000,   # loci used to estimate migration
               N_sim_mig         = 20,
               chains_mig        = 2,
               iter_warmup_mig   = 200,
               iter_sampling_mig = 200,
               iter_warmup_gea   = 100,
               iter_sampling_gea = 100,
               par = 12)                   # cores to use
str(res, max.level = 1)
```

Env/locus pair -log(e) p-values are located in `res$SNP_pvalues`. `res`
also contains a suite of plots and other outputs. See the vignette for
details.

# Large datasets

Whole-genome data sets are usually too large to fit in a single model
run. The recommended approach is to split the data into chunks of loci
and run each chunk separately, ideally in parallel on a computing
cluster.

mechGEA has a few utilities to help with this:

1.  [`get_chunks()`](https://hemstrow.github.io/mechGEA/reference/get_chunks.md)
    splits a “.beagle” or “.beagle.gz” file into chunks containing n
    loci.
2.  `example_slurm_job` generates directories, data and metadata files,
    and a ready-to-run job script for SLURM clusters (which can be
    modified for other job-management engines). This script will run
    `mechGEA` on each BEAGLE chunk.

## Splitting data into chunks

[`get_chunks()`](https://hemstrow.github.io/mechGEA/reference/get_chunks.md)
splits a BEAGLE file into smaller BEAGLE files of `loci_per_chunk` loci
each, reading the input in pieces so that the full file never has to be
loaded into memory. e.g recommended 50k SNPs per chunk at max for whole
genome dataset

``` r

bfile  <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
chunks <- get_chunks(bfile, 
                     loci_per_chunk = 500, # loci per chunk, use 50-100k for WGS
                     outdir = ".") # save files to pwd. Default uses beagle dir.
chunks
#> [1] "chunk_0.beagle.gz" "chunk_1.beagle.gz" "chunk_2.beagle.gz"
```

## Running on a SLURM cluster

[`example_slurm_job()`](https://hemstrow.github.io/mechGEA/reference/example_slurm_job.md)
generates a SLURM array job that runs
[`mechGEA()`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md)
on each BEAGLE file as a separate task:

``` r

sh <- example_slurm_job(beagle_files    = chunks,
                        sample_metadata = sample_metadata,
                        output_dir      = "results",
                        conda_env       = "mechGEA")
sh
```

[`example_slurm_job()`](https://hemstrow.github.io/mechGEA/reference/example_slurm_job.md)
takes all of the arguments that `mechGEA` uses, as well as resource and
env pruning arguments.

This writes the following files to `./control/` (change this with
`control_file_prefix`):

| File | Purpose |
|----|----|
| `example_slurm_array_job.sh` | SLURM array script, one task per BEAGLE file |
| `example_array_rscript.R` | R script run by each task |
| `beaglelist.txt` | list of BEAGLE files |
| `sample_metadata.RDS` | prepared metadata (`patch_info`, `patches_xy`, `env`) |
| `parmfile.RDS` | additional arguments passed to [`mechGEA()`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md) |

Each task converts its BEAGLE file with
[`beagle2pl()`](https://hemstrow.github.io/mechGEA/reference/beagle2pl.md),
runs
[`mechGEA()`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md)
in its own sub-directory of `output_dir`
(`partition_<task>_<chunk name>`), and writes its results there. The
final outputs (e.g. the results of `res <- mechGEA(...)`) will be saved
as `final_results.RDS`.

To run, edit the `example_slurm_array_job.sh` for your specific needs
and execute.

``` bash
sbatch example_slurm_array_job.sh
```

## Citation

If you use mechGEA, please also cite:

> Goel, N., Bossu, C.M., Van Ee, J.J., Zavaleta, E., Ruegg, K.C. &
> Hooten, M.B. (2025). Identifying genomic adaptation to local climate
> using a mechanistic evolutionary model. *Methods in Ecology and
> Evolution*, 16(10), 2448–2460.
> <https://doi.org/10.1111/2041-210X.70117>

> Goel, N., Bossu, C.M., Yi, S., Robertson, E.C., Brown, T.M., Bolton,
> P.E., Vernasco, B.J., Amirkhiz, R.G., Zavaleta, E., Ruegg, K.C. &
> Hooten, M.B. (2026). Identifying adaptive variation in spatially
> structured populations using low-coverage whole-genome sequencing
> data. *Molecular Biology and Evolution*, 43(10), msag255.
> <https://doi.org/10.1093/molbev/msag255>
