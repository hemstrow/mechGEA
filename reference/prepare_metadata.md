# Format metadata for analysis.

From a [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object
containing patch and individual IDs as well as environmental variables,
create the three metadata input files required by
[`mechGEA`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md). An
example input is in
[`sample_metadata`](https://hemstrow.github.io/mechGEA/reference/sample_metadata.md).
An example output is in
[`example_gl`](https://hemstrow.github.io/mechGEA/reference/example_gl.md).

## Usage

``` r
prepare_metadata(
  sample_meta,
  env_mean_method = "unweighted",
  idp = 1,
  resort = FALSE
)
```

## Arguments

- sample_meta:

  [`sf`](https://r-spatial.github.io/sf/reference/sf.html) POINT
  geometry containing two required columns: `patch` containing patch IDs
  and `ind_id` containing individual IDs. Additional columns specify
  named environmental variables. Patches must be in continuous blocks
  and sorted identically to the likelihood data also provided to
  [`mechGEA`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md).
  NOTE: `.beagle` or `.beagle.gz` files can be re-sorted on import using
  the
  [`beagle2pl`](https://hemstrow.github.io/mechGEA/reference/beagle2pl.md)
  function.

- env_mean_method:

  character, default "unweighted". Specifies how environmental variables
  are averaged across samples. Options:

  - "unweighted": direct means across samples within patches.

  - "weighted": inverse-distance scaled averages. Scaling power given by
    `idp`.

- idp:

  integer, default 1. Inverse-distance power used for weighted
  environmental averages. Ignored if `env_mean_method == "weighted"`.
  Higher values give less weight to points distant from the patch
  centroid.

- resort:

  logical, default FALSE. If TRUE, outputs will be automatically
  resorted such that patches are in blocks if needed. An additional
  output (n_ord) will be returned with the order used to resort that can
  be passed to
  [`beagle2pl`](https://hemstrow.github.io/mechGEA/reference/beagle2pl.md)
  to resort genotypes identically on read. WARNING: failing to resort
  beagle on read will cause incorrect downstream results.

## Author

William Hemstrom
