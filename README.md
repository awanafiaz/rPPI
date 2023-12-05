
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rPPI

<!-- badges: start -->

[![R-CMD-check](https://github.com/awanafiaz/rPPI/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/awanafiaz/rPPI/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## <img src="man/figures/PPI_logo.png" align="right" height="200" style="float:right; height:200px;">

Increased emphasis on reproducibility in study design and statistical
analysis has both prompted and necessitated national efforts to make
large-scale data sources widely available. With nationally
representative survey data available now more than ever, researchers can
answer more complex scientific questions with improved generalizability
and statistical precision. With the rise in interest in causal
inference, a branch of statistics dedicated to understanding dynamic
cause-and-effect mechanisms, there has been a corresponding shift in
focus toward making causal statements using complex survey data.

Drawing causal conclusions from observational data requires careful
considerations of both the appropriate quantity to estimate and the
rigor of statistical methods to achieve this end. The `svyate` package
introduces a novel strategy for incorporating sampling mechanisms such
as survey weights in a statistically valid analysis to make causal
statements when an individual’s probability of selection depends on the
treatment or characteristic of interest.

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
