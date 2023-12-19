
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BirdFlowExtras

<!-- badges: start -->
<!-- badges: end -->

**BirdFlowExtras** extends **BirdFlowR** with functions that are
anticipated to be used by fewer users than the core functions; or that
rely on packages that are not on CRAN or Bioconductor and thus cannot be
used by packages on CRAN.

## Installation

You can install the development version of BirdFlowExtras like so:

``` r
if(!require("remotes"))
  install.packages("remotes") 
remotes::install_github("birdflow-science/BirdFlowR")
```

## Example

Calculates Migratory Connectivity (MC) for a BirdFlowR model. This
relies on the **MigConnectivity** package to do most of the work.

``` r
library(BirdFlowExtras)
bf <- BirdFlowModels::amewoo
calc_birdflow_mc(bf, season = "postbreeding")
#> [1] 0.2690275
```
