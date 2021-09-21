
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SpaceX

<!-- badges: start -->
<!-- badges: end -->

The goal of SpaceX is to provide shared and cluster specfic gene
co-expression networks for spatial transcriptomics data.

## Installation

You can install the released version of SpaceX from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("SpaceX")
```

## Example

This is a basic example which shows you how to solve a common problem:

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
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
