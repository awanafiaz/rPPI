
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rPPI

<!-- badges: start -->

[![R-CMD-check](https://github.com/awanafiaz/rPPI/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/awanafiaz/rPPI/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## <img src="man/figures/PPI_logo.png" align="right" height="200" style="float:right; height:200px;">

With the rapid advancement of artificial intelligence (AI) and machine
learning (ML) methods and owing to financial and domain-specific
restrictions, researchers from wide variety of disciplines are now
increasingly using predictions from pre-trained algorithms as outcome
variables in statistical analyses. However, such practices have been
shown to produce biased estimates and misleading inference [Wang et al.,
2020](https://www.pnas.org/doi/suppl/10.1073/pnas.2001238117). The
statistical challenges encountered in post-prediction inference (PPI)
include (i) correct specification of the relationship between predicted
outcomes and their true unobserved counterparts, (ii) robustness of the
ML models to resampling or uncertainty about the training data, and
(iii) appropriately propagating both bias and uncertainty from
predictions into downstream inferential tasks.

Since the seminal work from [Wang et
al.](https://www.pnas.org/doi/suppl/10.1073/pnas.2001238117), published
in only 2020, multiple newer methods have been developed in quick
succession as researchers are increasingly putting more focus on this
new paradigm of statistical inference. Since each methods has their own
set of advantages and restrictions, we develop the `rPPI` package to
provide logistical software support to practitioners by allowing them to
utilize all extant methods under the umbrella of a single `r`package.

To make the utilization of the package convenient for users, we provide
guidance on installation and use of the package and its functions in the
following:

## Installation

You can install the development version of rPPI from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("awanafiaz/rPPI")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(rPPI)
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
up-to-date. `devtools::build_readme()` is handy for this.

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
