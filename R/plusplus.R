#===============================================================================
#
#  PROGRAM: plusplus.R (PPI++)
#
#  AUTHORS: Kentaro Hoffman (khoffm3@uw.edu)
#           Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: Implementation of various algorithms from Angelopoulos et al. (2023)
#
#           PPI++
#
#  INPUTS:  N/A
#
#  OUTPUS:  Prediction-powered inference functions for various target estimands: # SS: Wrapper function around these?
#
#           1. ppi_ols_ci: Linear Regression                                     # SS: Better naming conventions?
#
#  Notes:   1. How do we want to handle the data argument? One argument with a
#              stacked dataset and a column that gives the set (tr, te, val)?
#              Or two dataset arguments, an unlabeled and labeled?
#
#  Updated: 2023-12-02
#
#===============================================================================

#=== HELPER FUNCTIONS FOR PPI++ ================================================

#--- NORMAL CONFIDENCE INTERVALS -----------------------------------------------

zconfint_generic <- function(mean, std_mean, alpha, alternative) {

  if (alternative %in% c("two-sided", "2-sided", "2s")) {

    zcrit <- qnorm(1 - alpha / 2)
    lower <- mean - zcrit * std_mean
    upper <- mean + zcrit * std_mean

  } else if (alternative %in% c("larger", "l")) {

    zcrit <- qnorm(alpha)
    lower <- mean + zcrit * std_mean
    upper <- Inf

  } else if (alternative %in% c("smaller", "s")) {

    zcrit <- qnorm(1 - alpha)
    lower <- -Inf
    upper <- mean + zcrit * std_mean

  } else {

    stop("Invalid alternative")
  }

  return(c(lower, upper))
}

#--- ESTIMATE POWER TUNING PARAMETER -------------------------------------------

calc_lhat_glm <- function(grads, grads_hat, grads_hat_unlabeled, inv_hessian,

  coord = NULL, clip = F) {

  n <- nrow(grads)
  N <- nrow(grads_hat_unlabeled)
  d <- ncol(inv_hessian)

  cov_grads <- matrix(0, nrow = d, ncol = d)

  for (i in 1:n) {

    cov_grads <- cov_grads + (1 / n) * (

      outer(grads[i,] - colMeans(grads), grads_hat[i,] - colMeans(grads_hat)) +

      outer(grads_hat[i,] - colMeans(grads_hat), grads[i, ] - colMeans(grads)))
  }

  var_grads_hat <- cov(rbind(grads_hat, grads_hat_unlabeled))

  if (is.null(coord)) {

    vhat <- inv_hessian

  } else {

    vhat <- inv_hessian %*% diag(d)[, coord]
  }

  if (d > 1) {

    num <- ifelse(is.null(coord),

      sum(diag(vhat %*% cov_grads %*% vhat)), vhat %*% cov_grads %*% vhat)

    denom <- ifelse(is.null(coord),

      2*(1 + n/N) * sum(diag(vhat %*% var_grads_hat %*% vhat)),

      2*(1 + n/N) * vhat %*% var_grads_hat %*% vhat)

  } else {

    num <- vhat * cov_grads * vhat

    denom <- 2*(1 + n/N) * vhat * var_grads_hat * vhat
  }

  lhat <- num / denom

  if (clip) {

    lhat <- pmax(0, pmin(lhat, 1))
  }

  return(as.numeric(lhat))
}

#--- ORDINARY LEAST SQUARES ----------------------------------------------------

ols <- function(X, Y, return_se = F) {

  regression <- lm(Y ~ X - 1)

  theta <- coef(regression)

  if (return_se) {

    se <- sqrt(diag(vcov(regression)))

    return(list(theta = theta, se = se))

  } else {

    return(theta)
  }
}

#--- WEIGHTED LEAST SQUARES ----------------------------------------------------

wls <- function(X, Y, w = NULL, return_se = F) {

  if (is.null(w) || all(w == 1)) {

    return(ols(X, Y, return_se = return_se))
  }

  regression <- lm(Y ~ X - 1, weights = w)

  theta <- coef(regression)

  if (return_se) {

    se <- sqrt(diag(vcov(regression)))

    return(list(theta = theta, se = se))

  } else {

    return(theta)
  }
}

#--- OLS GRADIENT AND HESSIAN --------------------------------------------------

ols_get_stats <- function(pointest, X, Y, Yhat, X_unlabeled, Yhat_unlabeled,

  w = NULL, w_unlabeled = NULL, use_unlabeled = T) {

  n <- nrow(X)
  N <- nrow(X_unlabeled)
  d <- ncol(X)

  if (is.null(w)) {

    w <- rep(1, n)

  } else {

    w <- w / sum(w) * n
  }

  if (is.null(w_unlabeled)) {

    w_unlabeled <- rep(1, N)

  } else {

    w_unlabeled <- w_unlabeled / sum(w_unlabeled) * N
  }

  hessian <- matrix(0, nrow = d, ncol = d)

  grads_hat_unlabeled <- matrix(0, nrow = N, ncol = ncol(X_unlabeled))

  if (use_unlabeled) {

    for (i in 1:N) {

      hessian <- hessian + w_unlabeled[i] / (N + n) *

        tcrossprod(X_unlabeled[i, ])

      grads_hat_unlabeled[i, ] <- w_unlabeled[i] * X_unlabeled[i, ] *

        (sum(X_unlabeled[i, ] * pointest) - Yhat_unlabeled[i])
    }
  }

  grads <- matrix(0, nrow = n, ncol = ncol(X))

  grads_hat <- matrix(0, nrow = n, ncol = ncol(X))

  for (i in 1:n) {

    if (use_unlabeled) {

      hessian <- hessian + w[i] / (N + n) * tcrossprod(X[i, ])

    } else {

      hessian <- hessian + w[i] / n * tcrossprod(X[i, ])
    }

    grads[i, ] <- w[i] * X[i, ] * (sum(X[i, ] * pointest) - Y[i])

    grads_hat[i, ] <- w[i] * X[i, ] * (sum(X[i, ] * pointest) - Yhat[i])
  }

  inv_hessian <- solve(hessian)

  return(list(grads = grads, grads_hat = grads_hat,

    grads_hat_unlabeled = grads_hat_unlabeled, inv_hessian = inv_hessian))
}

#--- POINT ESTIMATE ------------------------------------------------------------

ppi_ols_pointestimate <- function(X, Y, Yhat, X_unlabeled, Yhat_unlabeled,

  lhat = NULL, coord = NULL, w = NULL, w_unlabeled = NULL) {

  n <- nrow(X)
  d <- ncol(X)
  N <- nrow(X_unlabeled)

  if (is.null(w)) {

    w <- rep(1, n)

  } else {

    w <- w / sum(w) * n
  }

  if (is.null(w_unlabeled)) {

    w_unlabeled <- rep(1, N)

  } else {

    w_unlabeled <- w_unlabeled / sum(w_unlabeled) * N
  }

  use_unlabeled <- is.null(lhat) || lhat != 0

  if (is.null(lhat)) {

    imputed_theta <- wls(X_unlabeled, Yhat_unlabeled, w = w_unlabeled)

  } else {

    imputed_theta <- wls(X_unlabeled, lhat * Yhat_unlabeled, w = w_unlabeled)
  }

  if (is.null(lhat)) {

    rectifier <- wls(X, Y - Yhat, w = w)

  } else {

    rectifier <- wls(X, Y - lhat * Yhat, w = w)
  }

  ppi_pointest <- imputed_theta + rectifier

  if (is.null(lhat)) {

    stats <- ols_get_stats(ppi_pointest, as.matrix(X), Y, Yhat,

      as.matrix(X_unlabeled), Yhat_unlabeled, w = w, w_unlabeled = w_unlabeled,

      use_unlabeled = use_unlabeled)

    lhat <- calc_lhat_glm(stats$grads, stats$grads_hat,

      stats$grads_hat_unlabeled, stats$inv_hessian, coord, clip = T)

    return(ppi_ols_pointestimate(X, Y, Yhat, X_unlabeled, Yhat_unlabeled,

      lhat = lhat, coord = coord, w = w, w_unlabeled = w_unlabeled))

  } else {

    return(ppi_pointest)
  }
}

#=== LINEAR REGRESSION =========================================================

#' PPI++ Linear Regression using Angelopoulos et al. (2023) Analytic Form
#'
#' @description
#' A short description...
#'
#' @details
#' Additional details...
#'
#' @param rec_form A formula defining the rectifier model. This should be of
#' the form Y - Yhat ~ X, where Y is the name of the column corresponding to
#' the observed outcome in the labeled data, Yhat is the name of the column
#' corresponding to the predicted outcome in the labeled data, and X generally
#' corresponds to the features of interest (e.g., X1 + X2).
#'
#' @param inf_form A formula defining the inference model. This should be of
#' the form Yhat ~ X, where Yhat is the name of the column corresponding to the
#' predicted outcome in the unlabeled data, and X generally corresponds to the
#' features of interest (e.g., X1 + X2).
#'
#' @param dat data in the form of the simdat function output
#'
#' @param alpha scalar type I error rate for hypothesis testing - values in (0, 1); defaults to 0.05
#'
#' @param alternative alternative hypothesis for hypothesis testing - options include "one-sided" or "two-sided"; defaults to two-sided
#'
#' @param lhat scalar power tuning parameter - values in (0, 1); defaults to NULL and is estimated
#'
#' @param coord vector coordinates for gradient descent; defaults to NULL and is estimated
#'
#' @param w vector weights for weighted glm in labeled dataset; defaults to NULL and is estimated
#'
#' @param w_unlabeled vector weights for weighted glm in unlabeled dataset; defaults to NULL and is estimated
#'
#' @returns A list of outputs: estimate of inference model parameters and corresponding standard errors
#'
#' @examples
#'
#' rec_form <- Y - Yhat ~ X1
#'
#' inf_form <- Yhat ~ X1
#'
#' dat <- simdat()
#'
#' ppi_ols_ci(rec_form, inf_form, dat = dat)
#'
#' @import stats
#'
#' @export

ppi_ols_ci <- function(rec_form, inf_form, dat, alpha = 0.05,

  alternative = "two-sided", lhat = NULL, coord = NULL, w = NULL,

  w_unlabeled = NULL) {

  X <- model.matrix(rec_form, data = dat[dat$set == "tst",])

  Y <- matrix(dat[dat$set == "tst", all.vars(rec_form)[1]], ncol = 1)

  Yhat <- dat[dat$set == "tst", all.vars(rec_form)[2]]

  X_unlabeled <- model.matrix(inf_form, data = dat[dat$set == "val",])

  Yhat_unlabeled <- matrix(dat[dat$set == "val", all.vars(inf_form)[1]], ncol = 1)

  n <- nrow(X)

  d <- ncol(X)

  N <- nrow(X_unlabeled)

  if (is.null(w)) {

    w <- rep(1, n)

  } else {

    w <- w / sum(w) * n
  }

  if (is.null(w_unlabeled)) {

    w_unlabeled <- rep(1, N)

  } else {

    w_unlabeled <- w_unlabeled / sum(w_unlabeled) * N
  }

  use_unlabeled <- is.null(lhat) || lhat != 0

  ppi_pointest <- ppi_ols_pointestimate(X, Y, Yhat, X_unlabeled, Yhat_unlabeled,

    lhat = lhat, coord = coord, w = w, w_unlabeled = w_unlabeled)

  stats <- ols_get_stats(ppi_pointest, as.matrix(X), Y, Yhat,

    as.matrix(X_unlabeled), Yhat_unlabeled, w = w, w_unlabeled = w_unlabeled,

    use_unlabeled = use_unlabeled)

  if (is.null(lhat)) {

    lhat <- calc_lhat_glm(stats$grads, stats$grads_hat,

      stats$grads_hat_unlabeled, stats$inv_hessian, coord, clip = T)

    return(ppi_ols_ci(rec_form, inf_form, dat, alpha = alpha,

      alternative = alternative, lhat = lhat, coord = coord, w = w,

      w_unlabeled = w_unlabeled))
  }

  var_unlabeled <- cov(lhat * stats$grads_hat_unlabeled)

  var <- cov(stats$grads - lhat * stats$grads_hat)

  Sigma_hat <- stats$inv_hessian %*%

    (n/N * var_unlabeled + var) %*% stats$inv_hessian

  # return(zconfint_generic(ppi_pointest, sqrt(diag(Sigma_hat) / n),
  #
  #   alpha = alpha, alternative = alternative))

  return(list(est = ppi_pointest, se = sqrt(diag(Sigma_hat) / n)))
}

#=== END =======================================================================
