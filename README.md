
<style>
body {
text-align: justify}
</style>
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
include: (1) correct specification of the relationship between predicted
outcomes and their true unobserved counterparts, (2) robustness of the
ML models to resampling or uncertainty about the training data, and (3)
appropriately propagating both bias and uncertainty from predictions
into downstream inferential tasks.

Since the seminal work from [Wang et
al.](https://www.pnas.org/doi/suppl/10.1073/pnas.2001238117), published
in only 2020, multiple newer methods have been developed in quick
succession as researchers are increasingly putting more focus on this
new paradigm of statistical inference. Since each methods has their own
set of advantages and restrictions, we have developed the `rPPI` package
to provide software support to researchers by allowing them to utilize
all extant methods under the umbrella of a single r-package.

To make the utilization of the package convenient for users, we provide
guidance on installation and use of the package and its functions in the
following:

## Installation

You can install the development version of rPPI from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")   ## If devtools is not already installed
devtools::install_github("awanafiaz/rPPI")
```

## Usage

We provide a simple example to demonstrate the basic use of the
functions included in `rPPI`. We frame our example in the following
manner to build an unifying example to be used for all the methods.

(i). Assume that we have access to a well-performing accurate AI/ML/DL
algorithm $f_{\text{x}}(\cdot)$ that can predict our outcome of interest
$Y$.

(ii). Next, consider that we have 2 data sets, a labeled data set, which
we call the **test set** $(X_{te}, Y_{te})$ and an unlabeled data set,
which we call the **validation set** $(X_{val)}$. Typically the the
labeled test set will be smaller in size compared to the unlabeled
validation set. Here we will consider them to be equal for brevity.

- We consider $X = (X_1, X_2, X_3, X_4)$ and $y$ is a scalar.
- The true data generating mechanism is
  $Y = \beta_1X_1 + \beta_2 X_2 + \beta_3 \ g(X_3) + \beta_4 \ g(X_4) + \epsilon,$
  where, $\epsilon = N(0, 1)$ and $g(\cdot)$ refers to some smoothing
  function.
- We specify,
  $\beta = \begin{bmatrix} \beta_1 \\ \beta_2 \\ \beta_3 \\ \beta_4 \end{bmatrix} = \begin{bmatrix} 1 \\ 0.5 \\ 3 \\ 4 \end{bmatrix}$.

(iii). Our interest is in performing inference on $H_0: \beta_1 = 0$ vs
$H_1: \beta_1 \ne 0$. That is our inference model is
$Y_{val} = \beta_0^* + \beta_1^* X_{val} + \epsilon^*$.

(iv). However, we do not observe $Y_{val}$. We instead only have access
to the predicted $\hat Y_{val} = f(X_{val})$.

We will now obtain the estimates from each different method.

``` r
## Load the library
library(rPPI)

## generate the data
set.seed(2023)
dat <- simdat(n = c(300, 300, 300), beta1 = 1)
```

#### 1. Original analytic method from Wang et al. (2020)

``` r
# Requires the specification of 
## 1. relationship model between observed y and predicted Y-hat 
rel_form <- Y ~ Yhat  ## we consider a basic linear function
## 2. inference model
inf_form <- Yhat ~ X1

rPPI::wang_analytic(rel_form, inf_form, dat = dat)
#> $est
#> [1] 11.561038  0.955248
#> 
#> $se
#> [1] 0.2401493 0.1706713
```

#### 2. Original bootstrap method from Wang et al. (2020)

``` r
# Requires the specification of 
## 1. relationship model between observed y and predicted Y-hat 
rel_form <- Y ~ Yhat  ## we consider a basic linear function
## 2. inference model
inf_form <- Yhat ~ X1
## 3. we also need to specify the number of bootstraps 
nboot <- 200
rPPI::wang_boot(rel_form, inf_form, dat = dat, nboot)
#> $est
#> [1] 11.5462734  0.9444002
#> 
#> $se_par
#> [1] 0.2394565 0.1706623
#> 
#> $se_npar
#> [1] 0.2364742 0.1807781
# the function returns both parametric (par) and non-parametric (npar) estimate of std.error (se)
```

#### 3. Prediction-powered inference method (Angelopoulos et al., 2023)

``` r
# Requires the specification of 
## 1. rectifier model
rec_form <- Y - Yhat ~ X1 
## 2. inference model
inf_form <- Yhat ~ X1

rPPI::angelopoulos_original(rec_form, inf_form, dat)
#> $est
#> [1] 11.697627  0.817898
#> 
#> $se
#> [1] 0.2357630 0.1790149
```

#### 4. Multiple-imputation method from Leek et al., (2023)

``` r
# Requires the specification of 
## 1. relationship model between observed y and predicted Y-hat 
rel_form <- Y ~ Yhat  ## we consider a basic linear function
## 2. inference model
inf_form <- Yhat ~ X1
m <- 100

rPPI::wang_mi(rel_form, inf_form, dat = dat, m)
#> $est
#> [1] 11.5855979  0.9351817
#> 
#> $se
#> [1] 0.3228812 0.2336300
```

## Vignette

For more advanced users and researchers, we provide more use cases and
examples in the package `vignettes`.

``` r
vignette("rPPI")
```

This will provide an extensive tutorial on `rPPI` and will demonstrate
how to use it.

## Feedback

For questions and comments or any other feedback, please contact the
developers.
