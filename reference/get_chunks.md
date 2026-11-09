# Chunk a beagle file into partitions for processing.

Creates n chunks of a beagle file, where sequential lines are sent to
different files round-robin style. Note that this is much faster on
Linux than on Windows, since it avoids reading anything into memory
using unix commands.

## Usage

``` r
get_chunks(
  beagle_file,
  loci_per_chunk = 10000,
  max_read = 50000,
  outdir = dirname(beagle_file)
)
```

## Arguments

- beagle_file:

  path to `.beagle` or `.beagle.gz` file.

- loci_per_chunk:

  numeric, default 10000. Number of loci to distribute to each chunk.

- max_read:

  numeric, default 50000. Number of lines to read for splitting at once.
  Ignored on Linux.

- outdir:

  path to output directory, default `dirname(beagle_file)` (the
  directory containing the input file).

## Value

A vector of filenames for the chunks.

## Author

William Hemstrom

## Examples

``` r
# generate 100 SNP chunks
beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
chunks <- get_chunks(beagle, loci_per_chunk = 100)
#> Splitting into 15 files with approx. 100 snps per file.
file.remove(chunks) # remove the chunked files for the example only.
#>  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
```
