
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SpaceX

<!-- badges: start -->
<!-- badges: end -->

The goal of SpaceX is to provide shared and cluster specfic gene
co-expression networks for spatial transcriptomics data.

## Installation

You can install the released version of SpaceX from
(<https://github.com/SatwikAch/SpaceX>) with:

``` r
devtools::install_github("SatwikAch/SpaceX")
```

``` r
library(SpaceX)
#> Loading required package: PQLseq
#> Registered S3 methods overwritten by 'robust':
#>   method              from      
#>   plot.covfm          fit.models
#>   print.covfm         fit.models
#>   summary.covfm       fit.models
#>   print.summary.covfm fit.models
#> rlm is already registered in the fit.models registry
#> covfm is already registered in the fit.models registry
```

``` r
## Reading the Breast cancer data

## Spatial locations
head(BC_loc)

## Gene expression for data
head(BC_count) 

## Data processing
G <-dim(BC_count)[2] ## number of genes
N <-dim(BC_count)[1] ## number of locations

## Application to SpaceX algorithm
BC_fit <- SpaceX(BC_count,BC_loc[,1:2],BC_loc[,3])

##Output
## SigmaPhi :: Shared Covariance matrix
## SigmaLambda :: Cluster specific Covaraince matrices
```

You can view the supplementary file at this link: https://bookdown.org/satwik91/SpaceX_supplementary/.
