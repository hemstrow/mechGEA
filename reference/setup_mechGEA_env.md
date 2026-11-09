# Set-up a conda environment with all `mechGEA` dependencies properly installed.

[`mechGEA`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md) can
be difficult to properly install given very exact install order
requirements to prevent conflicts and bugs on some systems. This
function sets up a properly-installed conda/mamba environment for future
use.

## Usage

``` r
setup_mechGEA_env(name = "mechGEA", conda_path = "DEFAULT")
```

## Arguments

- name:

  character, default "mechGEA". Name of the environment to create.

- conda_path:

  character, default "DEFAULT". Path to conda binary. If "DEFAULT",
  [`conda_binary`](https://rstudio.github.io/reticulate/reference/conda-tools.html)
  will be used to attempt a guess and install if no conda install
  located.

## Value

Generates a conda environment with
[`mechGEA`](https://hemstrow.github.io/mechGEA/reference/mechGEA.md) and
all dependencies correctly installed. Load with `conda activate mechGEA`
(or whatever name was used during setup) called from the command line.

## Author

Eric Robertson

William Hemstrom

## Examples

``` r
# setup an environment using the default conda path
if (FALSE) { # \dontrun{
setup_mechGEA_env("mechGEA")

# then, in command line (not R), run:
# conda activate mechGEA
# Rscript /path/to/my/mechGEA/analysis/Rscript.R
} # }
```
