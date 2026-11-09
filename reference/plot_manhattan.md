# Generate a manhattan plot one or more mechGEA outputs.

Creates a ggplot-based manhattan plot, where chromosomes/scaffolds/etc
are concatenated along the x-axis. Works with
[`mechGEA`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md)
SNP_pvalue outputs; take either a single output or a list of outputs
from multiple chunks (see details).

## Usage

``` r
plot_manhattan(
  SNP_pvalues,
  SNP_metadata,
  vlines = FALSE,
  vline_width = 0.25,
  significant = 5,
  suggestive = 3,
  colors = c("black", "slategray3"),
  chr_order = NULL,
  abbreviate_labels = FALSE,
  env_colors = NULL
)
```

## Arguments

- SNP_pvalues:

  data.frame or similar like that produced by
  [`mechGEA`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md)
  in `$SNP_pvalues` OR the output from
  [`mechGEA`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md)
  \\
  [`run_stan_GEA_model`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md)
  OR a list of a either of the two (mixing OK). If the former, contains
  three columns named "SNP", "Env", and "nlogp" containing SNP IDs,
  names for environmental effects, and -log p-values for associations.

- SNP_metadata:

  data.frame or similar containing metadata on SNPs OR a path to a
  ".beagle" or ".beagle.gz" file containing SNP metadata OR a list of a
  either of the two (mixing OK). If the former, must have three columns
  named "SNP", "CHROM", and "POS" containing SNP IDs, chromosomes, and
  positions for each SNP. Each SNP ID in the `SNP_pvalues` "SNP" column
  must be represented in the "SNP" column of `SNP_metadata`. Assumes
  sorting is in the same order as in `SNP_pvalues`, across lists or not.

- vlines:

  character (color) or FALSE, default FALSE. If a color, vertical
  separator lines will be drawn between each chromosome. Widths
  controlled by `vline_width`.

- vline_width:

  numeric, default 2. Width of chromosome separator lines. Ignored if
  `vlines` is FALSE.

- significant:

  numeric, default NULL. Value at which a line will be drawn designating
  significant SNPs. If highlight = "significant", SNPs above this level
  will also be labeled.

- suggestive:

  numeric, default NULL. Value at which a line will be drawn designating
  suggestive SNPs. If highlight = "suggestive", SNPs above this level
  will also be labeled.

- colors:

  character, default c("black", "slategray3"). Colors to alternate
  across chromosomes.

- chr_order:

  character, default NULL. If provided, an ordered vector of
  chromosome/scaffold/etc names by which to sort output.

- abbreviate_labels:

  numeric or FALSE, default FALSE. If a numeric value, x-axis chromosome
  names will be abbreviated using
  [`abbreviate`](https://rdrr.io/r/base/abbreviate.html), with each
  abbreviated label having the minimum length specified. Helpful when
  chromosome/scaffold/etc names are very long.

- env_colors:

  character, default NULL. An optional vector of colors to use for each
  environmental variable.

## Author

William Hemstrom

## Examples

``` r
if (FALSE) { # \dontrun{
# run a model
beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
phred <- beagle2pl(beagle, lines = 1:150)
input_data <- prepare_metadata(sample_metadata)
res <- mechGEA(phred,
               patch_info = input_data$patch_info,
               patches_xy = input_data$patches_xy,
               migration_loci = 50,
               env = input_data$env,
               iter_sampling_gea = 100,
               iter_warmup_gea = 200,
               iter_sampling_mig = 100,
               iter_warmup_mig = 200)

# plot results
plot_manhattan(res, beagle)
## or
plot_manhattan(res$SNP_pvalues, beagle)
# or could provide a data.frame of SNP metadata


# combine across result runs. Note that using chunk_beagle() is easier, than
# manually splitting, but for the example we're just making two chunks
beagle <- system.file("extdata", "example.beagle.gz", package = "mechGEA")
system(paste0("zcat ", beagle2,
              " | head -n ", 301, " | gzip > test_sub.beagle.gz"))
b <- fread(" test_sub.beagle.gz")
b <- b[151:300,] # the second batch of lines
fwrite(b, "test_sub.beagle.gz", sep = "\t",
       col.names = TRUE, row.names = FALSE)

phred <- beagle2pl("test_sub.beagle.gz")
res2 <- mechGEA(phred,
                patch_info = input_data$patch_info,
                patches_xy = input_data$patches_xy,
                migration_loci = 50,
                env = input_data$env,
                iter_sampling_gea = 100,
                iter_warmup_gea = 200,
                iter_sampling_mig = 100,
                iter_warmup_mig = 200)

# plot results
plot_manhattan(list(res, res2), list(beagle, "test_sub.beagle.gz"))
# Any combination of res/res$SNP_pvalues and filepaths/tables of metadata
# are allowed in lists. Assumes everything IS IN THE SAME ORDER!
} # }
```
