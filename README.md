
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

## Here we are randomly choosing 20 genes from 40 locations to run our SpaceX algorithm.
set.seed(123)
samp_G <- sample(1:G,20,replace = F)
samp_N <- sample(1:N,40,replace = F)

BC_loc <- BC_loc[samp_N,]
BC_count <- BC_count[samp_N,samp_G]


## Application to SpaceX algorithm
BC_fit <- SpaceX(BC_count,BC_loc)

## SigmaPhi :: Shared Covariance matrix
## SigmaLambda :: Cluster specific Covaraince matrices
```
