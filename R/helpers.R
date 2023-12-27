#===============================================================================
# HELPER FUNCTIONS
#===============================================================================




#=== PPI++ =====================================================================

#--- ORDINARY LEAST SQUARES ----------------------------------------------------

#' Ordinary Least Squares
#'
#' @description
#' Computes the ordinary least squares coefficients.
#'
#' @param X (matrix): n x p matrix of covariates.
#'
#' @param Y (vector): p-vector of outcome values.
#'
#' @param return_se (bool, optional): Whether to return the standard errors of
#' the coefficients.
#'
#' @returns (list): A list containing the following:
#'
#' \describe{
#'    \item{theta}{(vector): p-vector of ordinary least squares estimates of
#'    the coefficients.}
#'    \item{se}{(vector): If return_se == TRUE, return the p-vector of
#'    standard errors of the coefficients.}
#' }
#'
#' @examples
#'
#' n <- 1000
#'
#' X <- rnorm(n, 1, 1)
#'
#' Y <- X + rnorm(n, 0, 1)
#'
#' ols(X, Y, return_se = T)
#'
#' @export

ols <- function(X, Y, return_se = F) {

  fit <- lm(Y ~ X - 1)

  theta <- coef(fit)

  if (return_se) {

    se <- sqrt(diag(vcov(fit)))

    return(list(theta = theta, se = se))

  } else {

    return(theta)
  }
}

#--- WEIGHTED LEAST SQUARES ----------------------------------------------------

#' Weighted Least Squares
#'
#' @description
#' Computes the weighted least squares estimate of the coefficients.
#'
#' @param X (matrix): n x p matrix of covariates.
#'
#' @param Y (vector): p-vector of outcome values.
#'
#' @param w (vector, optional): n-vector of sample weights.
#'
#' @param return_se (bool, optional): Whether to return the standard errors of
#' the coefficients.
#'
#' @returns (list): A list containing the following:
#'
#' \describe{
#'    \item{theta}{(vector): p-vector of weighted least squares estimates of
#'    the coefficients.}
#'    \item{se}{(vector): If return_se == TRUE, return the p-vector of
#'    standard errors of the coefficients.}
#' }
#'
#' @examples
#'
#' n <- 1000
#'
#' X <- rnorm(n, 1, 1)
#'
#' w <- rep(1, n)
#'
#' Y <- X + rnorm(n, 0, 1)
#'
#' wls(X, Y, w = w, return_se = T)
#'
#' @export

wls <- function(X, Y, w = NULL, return_se = F) {

  if (is.null(w) || all(w == 1)) {

    return(ols(X, Y, return_se = return_se))
  }

  fit <- lm(Y ~ X - 1, weights = w)

  theta <- coef(fit)

  if (return_se) {

    se <- sqrt(diag(vcov(fit)))

    return(list(theta = theta, se = se))

  } else {

    return(theta)
  }
}















################################################################################
################################################################################




#--- OLS GRADIENT AND HESSIAN --------------------------------------------------

"""Computes the statistics needed for the OLS-based prediction-powered inference.

    Args:
        pointest (ndarray): A point estimate of the coefficients.
        X (ndarray): Covariates for the labeled data set.
        Y (ndarray): Labels for the labeled data set.
        Yhat (ndarray): Predictions for the labeled data set.
        X_unlabeled (ndarray): Covariates for the unlabeled data set.
        Yhat_unlabeled (ndarray): Predictions for the unlabeled data set.
        w (ndarray, optional): Sample weights for the labeled data set.
        w_unlabeled (ndarray, optional): Sample weights for the unlabeled data set.
        use_unlabeled (bool, optional): Whether to use the unlabeled data set.

    Returns:
        grads (ndarray): Gradient of the loss function with respect to the coefficients.
        grads_hat (ndarray): Gradient of the loss function with respect to the coefficients, evaluated using the labeled predictions.
        grads_hat_unlabeled (ndarray): Gradient of the loss function with respect to the coefficients, evaluated using the unlabeled predictions.
        inv_hessian (ndarray): Inverse Hessian of the loss function with respect to the coefficients.
    """

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


#=== END =======================================================================











################################################################################
################################################################################

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

  if(is.null(dim(grads))) grads <- matrix(grads, ncol = 1)

  if(is.null(dim(grads_hat))) grads_hat <- matrix(grads_hat, ncol = 1)

  if(is.null(dim(grads_hat_unlabeled))) grads_hat_unlabeled <- matrix(grads_hat_unlabeled, ncol = 1)

  n <- nrow(grads)
  N <- nrow(grads_hat_unlabeled)
  d <- ncol(inv_hessian)

  cov_grads <- matrix(0, nrow = d, ncol = d)

  for (i in 1:n) {

      cov_grads <- cov_grads + (1 / n) * (

        outer(grads[i,] - colMeans(grads), grads_hat[i,] - colMeans(grads_hat)) +

          outer(grads_hat[i,] - colMeans(grads_hat), grads[i,] - colMeans(grads)))

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





