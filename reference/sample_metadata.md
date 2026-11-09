# Example sample metadata

A simulated dataset containing patch data, lat/lon coordinates, and
environmental data for the 200 individuals in
[`example_gl`](https://hemstrow.github.io/mechGEA/reference/example_gl.md),
in an [`sf`](https://r-spatial.github.io/sf/reference/sf.html) format as
expected by
[`prepare_metadata`](https://hemstrow.github.io/mechGEA/reference/prepare_metadata.md).

## Usage

``` r
sample_metadata
```

## Format

An object of class `sf` (inherits from `data.frame`) with 200 rows and 8
columns.

## Details

Individuals are in rows, sorted identically to genotypic data. Expects
columns named:

- patch: name of the patch from which an individual sample was
  collected.

- ind_id: individual IDs.

All remaining variables must be numeric and are assumed to be
uncorrelated environmental variables. The function
[`prune_envs`](https://hemstrow.github.io/mechGEA/reference/prune_envs.md)
can be used to generate uncorrelated variables from a full matrix of
envs.
