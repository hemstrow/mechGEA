# Check that the mechGEA package is fully, correctly installed.

This function runs a series of checks to see that mechGEA is correctly
installed. It runs automatically during the
[`mechGEA`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md)
pipeline, but can also be run manually. It only needs to be run once per
session.

## Usage

``` r
check_correct_install(status = Sys.getenv("mechGEA_full_install"))
```

## Arguments

- stats:

  install status, fetched from a system environment (which should not be
  manually updated; this function will check and update the status as
  needed).

## Author

William Hemstrom
