# Format a beagle file as phred-scaled likelihoods.

Takes a beagle formatted input file and converts it to a matrix of
phred- scaled likelihoods, as needed by `mechGEA`. See details.

## Usage

``` r
beagle2pl(file, lines = NULL, ind_order = NULL, ind_names = NULL)
```

## Arguments

- file:

  character. File path to the beagle or beagle.gz file containing the
  raw data.

- lines:

  numeric, default NULL. Optional vector of numbers noting which lines
  in the beagle file to keep. This must be sequential due to
  [`fread`](https://rdrr.io/pkg/data.table/man/fread.html)'s
  restrictions.

- ind_order:

  numeric, default NULL. Numeric vector indicating final sorted order to
  use for individuals the output matrix (`c(2,1,3)`, for example, would
  put the second individual in `file` first, then the first individual
  second, and then the third individual last in the final matrix, for
  example). If NULL, the order will be unchanged.

- ind_names:

  character, default NULL. Names for individuals in the final matrix. If
  `ind_order` is provided, this should match the initial order of sample
  names, not the final order.

## Value

A \\p \times n\\ matrix, where p is the number of loci and n is the
number of individuals. See details for format.

## Details

The input file format required by this package has one row for each
locus and one column for each individual. Phred-scaled likelihoods (PL)
for the three possible genotypes at each locus/individual are comma
seperated (major homozygote, heterozygote, minor heterozygote). A higher
PL constitutes a lower likelihood; by convention the minimum PL is
subtracted from each triplet such that the most likely genotype has a PL
of 0.

Phred-likelihoods can be calculated according to:
\$\$PL(\text{genotype}) = -10 \times \log\_{10}\left( P(\text{Data} \mid
\text{Genotype}) \right)\$\$

Where \\P(\text{Data} \mid \text{Genotype})\\ is the a genotype
likelihood as produced by something like `ANGSD`'s `doGlf`.

## Author

William Hemstrom

## Examples

``` r
# example code

beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
pl <- beagle2pl(beagle,
               lines = 1:100) # first 100 lines
```
