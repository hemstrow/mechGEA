# Prune correlated (environmental) variables.

Conducts pruning of envirnonmental values according to correlation by
removing the variable correlated with the most other variables
iteratively until no correlated pairs remain. Ties removed at random.

## Usage

``` r
prune_envs(
  envs,
  r2 = 0.7,
  scale = TRUE,
  seed = floor(stats::runif(1, 1, 10000))
)
```

## Arguments

- envs:

  numeric matrix/data.frame/data.table of

- r2:

  numeric, default 0.7. Threshold to use for prunning; correlations
  ([`cor`](https://rdrr.io/r/stats/cor.html)) values higher than `r2`
  will result in the removal of envs.

- scale:

  logical, default TRUE. If TRUE, environmental vairables will be scaled
  to mean 0, sd 1. This is recommended for subsequent model fitting by
  [`mechGEA`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md).

- seed:

  integer, default a random number between 1 and 10,000. Seed to use for
  breaking ties.

## Author

William Hemstrom
