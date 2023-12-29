#===============================================================================
# PPI++ LOGISTIC REGRESSION
#===============================================================================

log1pexp <- function(x) {

  idxs <- x > 10

  out <- numeric(length(x))

  out[idxs] <- x[idxs]

  out[!idxs] <- log1p(exp(x[!idxs]))

  return(out)
}

#=== PPI++ LOGISTIC REGRESSION POINT ESTIMATE ==================================

# """Computes the prediction-powered point estimate of the logistic regression coefficients.
#
# Args:
#   X (ndarray): Covariates corresponding to the gold-standard labels.
# Y (ndarray): Gold-standard labels.
# Yhat (ndarray): Predictions corresponding to the gold-standard labels.
# X_unlabeled (ndarray): Covariates corresponding to the unlabeled data.
# Yhat_unlabeled (ndarray): Predictions corresponding to the unlabeled data.
# lhat (float, optional): Power-tuning parameter (see `[ADZ23] <https://arxiv.org/abs/2311.01453>`__). The default value `None` will estimate the optimal value from data. Setting `lhat=1` recovers PPI with no power tuning, and setting `lhat=0` recovers the classical point estimate.
# coord (int, optional): Coordinate for which to optimize `lhat`. If `None`, it optimizes the total variance over all coordinates. Must be in {1, ..., d} where d is the shape of the estimand.
# optimizer_options (dict, optional): Options to pass to the optimizer. See scipy.optimize.minimize for details.
# w (ndarray, optional): Sample weights for the labeled data set.
# w_unlabeled (ndarray, optional): Sample weights for the unlabeled data set.
#
# Returns:
#   ndarray: Prediction-powered point estimate of the logistic regression coefficients.
#
# Notes:
#   `[ADZ23] <https://arxiv.org/abs/2311.01453>`__ A. N. Angelopoulos, J. C. Duchi, and T. Zrnic. PPI++: Efficient Prediction Powered Inference. arxiv:2311.01453, 2023.
# """

ppi_logistic_pointestimate <- function(X, Y, Yhat, X_unlabeled, Yhat_unlabeled,

  lhat = NULL, coord = NULL, optimizer_options = NULL,

  w = NULL, w_unlabeled = NULL) {

  n <- nrow(Y)
  d <- ncol(X)
  N <- nrow(Yhat_unlabeled)

  if (is.null(w)) w <- rep(1, n) else w <- w / sum(w) * n

  if (is.null(w_unlabeled)) w_unlabeled <- rep(1, N) else w_unlabeled <- w_unlabeled / sum(w_unlabeled) * N

  if (is.null(optimizer_options) || !("ftol" %in% names(optimizer_options))) {

    optimizer_options <- list(ftol = 1e-15)
  }

  theta <- coef(glm(Y ~ . - 1, data = data.frame(Y, X), family = binomial))

  theta <- matrix(theta, ncol = 1)

  lhat_curr <- ifelse(is.null(lhat), 1, lhat)

  rectified_logistic_loss <- function(theta) {

    sum(w_unlabeled * (-Yhat_unlabeled * (X_unlabeled %*% theta) +

      log1pexp(X_unlabeled %*% theta))) * lhat_curr / N -

      sum(w * (-Yhat * (X %*% theta) + log1pexp(X %*% theta))) * lhat_curr / n +

      sum(w * (-Y * (X %*% theta) + log1pexp(X %*% theta))) / n
  }

  rectified_logistic_grad <- function(theta) {

    X_unlabeled_t <- t(X_unlabeled)

    X_t <- t(X)

    lhat_curr / N * X_unlabeled_t %*% (w_unlabeled *

      (plogis(X_unlabeled %*% theta) - Yhat_unlabeled)) -

      lhat_curr / n * X_t %*% (w * (plogis(X %*% theta) - Yhat)) +

      1 / n * X_t %*% (w * (plogis(X %*% theta) - Y))
  }

  ppi_pointest <- optim(par = theta, fn = rectified_logistic_loss,

    gr = rectified_logistic_grad, method = "L-BFGS-B",

    control = list(ftol = optimizer_options$ftol))$par

  if (is.null(lhat)) {

    stats <- logistic_get_stats(ppi_pointest, X, Y, Yhat, X_unlabeled,

      Yhat_unlabeled, w, w_unlabeled)

    lhat <- calc_lhat_glm(stats$grads, stats$grads_hat,

      stats$grads_hat_unlabeled, stats$inv_hessian, clip = TRUE)

    return(ppi_logistic_pointestimate(X, Y, Yhat, X_unlabeled, Yhat_unlabeled,

      optimizer_options = optimizer_options, lhat = lhat, coord = coord, w = w,

      w_unlabeled = w_unlabeled))

  } else {
    return(ppi_pointest)
  }
}

#=== PPI++ LOGISTIC GRADIENT AND HESSIAN =======================================

# """Computes the statistics needed for the logistic regression confidence interval.
#
#     Args:
#         pointest (ndarray): Point estimate of the logistic regression coefficients.
#         X (ndarray): Covariates corresponding to the gold-standard labels.
#         Y (ndarray): Gold-standard labels.
#         Yhat (ndarray): Predictions corresponding to the gold-standard labels.
#         X_unlabeled (ndarray): Covariates corresponding to the unlabeled data.
#         Yhat_unlabeled (ndarray): Predictions corresponding to the unlabeled data.
#         w (ndarray, optional): Standard errors of the gold-standard labels.
#         w_unlabeled (ndarray, optional): Standard errors of the unlabeled data.
#         use_unlabeled (bool, optional): Whether to use the unlabeled data.
#
#     Returns:
#         grads (ndarray): Gradient of the loss function on the labeled data.
#         grads_hat (ndarray): Gradient of the loss function on the labeled predictions.
#         grads_hat_unlabeled (ndarray): Gradient of the loss function on the unlabeled predictions.
#         inv_hessian (ndarray): Inverse Hessian of the loss function on the unlabeled data.
#     """

logistic_get_stats <- function(pointest, X, Y, Yhat, X_unlabeled,

  Yhat_unlabeled, w = NULL, w_unlabeled = NULL, use_unlabeled = TRUE) {

  n <- nrow(Y)
  d <- ncol(X)
  N <- nrow(Yhat_unlabeled)

  if (is.null(w)) w <- rep(1, n) else w <- w / sum(w) * n

  if (is.null(w_unlabeled)) w_unlabeled <- rep(1, N) else w_unlabeled <- w_unlabeled / sum(w_unlabeled) * N

  mu <- plogis(X %*% pointest)

  mu_til <- plogis(X_unlabeled %*% pointest)

  hessian <- matrix(0, nrow = d, ncol = d)

  grads_hat_unlabeled <- matrix(0, nrow = N, ncol = d)

  if (use_unlabeled) {

    for (i in 1:N) {

      hessian <- hessian + w_unlabeled[i] / (N + n) * mu_til[i] * (1 - mu_til[i]) * tcrossprod(X_unlabeled[i, ])

      grads_hat_unlabeled[i, ] <- w_unlabeled[i] * X_unlabeled[i, ] * (mu_til[i] - Yhat_unlabeled[i])
    }
  }

  grads <- matrix(0, nrow = n, ncol = d)

  grads_hat <- matrix(0, nrow = n, ncol = d)

  for (i in 1:n) {

    if (use_unlabeled) {

      hessian <- hessian + w[i] / (N + n) * mu[i] * (1 - mu[i]) * tcrossprod(X[i, ])

    } else {

      hessian <- hessian + w[i] / n * mu[i] * (1 - mu[i]) * tcrossprod(X[i, ])
    }

    grads[i, ] <- w[i] * X[i, ] * (mu[i] - Y[i])

    grads_hat[i, ] <- w[i] * X[i, ] * (mu[i] - Yhat[i])
  }

  inv_hessian <- solve(hessian)

  return(

    list(grads = grads, grads_hat = grads_hat,

         grads_hat_unlabeled = grads_hat_unlabeled, inv_hessian = inv_hessian))
}

#=== PPI++ LOGISTIC REGRESSION =================================================

# """Computes the prediction-powered confidence interval for the logistic regression coefficients using the PPI++ algorithm from `[ADZ23] <https://arxiv.org/abs/2311.01453>`__.
#
#     Args:
#         X (ndarray): Covariates corresponding to the gold-standard labels.
#         Y (ndarray): Gold-standard labels.
#         Yhat (ndarray): Predictions corresponding to the gold-standard labels.
#         X_unlabeled (ndarray): Covariates corresponding to the unlabeled data.
#         Yhat_unlabeled (ndarray): Predictions corresponding to the unlabeled data.
#         alpha (float, optional): Error level; the confidence interval will target a coverage of 1 - alpha. Must be in the range (0, 1).
#         alternative (str, optional): Alternative hypothesis, either 'two-sided', 'larger' or 'smaller'.
#         lhat (float, optional): Power-tuning parameter (see `[ADZ23] <https://arxiv.org/abs/2311.01453>`__). The default value `None` will estimate the optimal value from data. Setting `lhat=1` recovers PPI with no power tuning, and setting `lhat=0` recovers the classical CLT interval.
#         coord (int, optional): Coordinate for which to optimize `lhat`. If `None`, it optimizes the total variance over all coordinates. Must be in {1, ..., d} where d is the shape of the estimand.
#         optimizer_options (dict, ooptional): Options to pass to the optimizer. See scipy.optimize.minimize for details.
#         w (ndarray, optional): Weights for the labeled data. If None, it is set to 1.
#         w_unlabeled (ndarray, optional): Weights for the unlabeled data. If None, it is set to 1.
#
#     Returns:
#         tuple: Lower and upper bounds of the prediction-powered confidence interval for the logistic regression coefficients.
#
#     Notes:
#         `[ADZ23] <https://arxiv.org/abs/2311.01453>`__ A. N. Angelopoulos, J. C. Duchi, and T. Zrnic. PPI++: Efficient Prediction Powered Inference. arxiv:2311.01453, 2023.
#     """

ppi_logistic_ci <- function(X, Y, Yhat, X_unlabeled, Yhat_unlabeled,

  alpha = 0.05, alternative = "two-sided", lhat = NULL, coord = NULL,

  optimizer_options = NULL, w = NULL, w_unlabeled = NULL) {

  n <- nrow(Y)
  d <- ncol(X)
  N <- nrow(Yhat_unlabeled)

  w <- if (is.null(w)) rep(1, n) else w / sum(w) * n

  w_unlabeled <- if (is.null(w_unlabeled)) rep(1, N) else w_unlabeled / sum(w_unlabeled) * N

  use_unlabeled <- is.null(lhat) || lhat != 0

  cat("use_unlabeled:", use_unlabeled, "\n")

  ppi_pointest <- ppi_logistic_pointestimate(X, Y, Yhat, X_unlabeled,

    Yhat_unlabeled, optimizer_options = optimizer_options, lhat = lhat,

    coord = coord, w = w, w_unlabeled = w_unlabeled)

  cat("ppi_pointest:", ppi_pointest, "\n")

  stats <- logistic_get_stats(ppi_pointest, X, Y, Yhat, X_unlabeled,

    Yhat_unlabeled, w, w_unlabeled, use_unlabeled = use_unlabeled)

  if (is.null(lhat)) {

    lhat <- calc_lhat_glm(stats$grads, stats$grads_hat,

      stats$grads_hat_unlabeled, stats$inv_hessian, clip = TRUE)

    cat("lhat:", lhat, "\n")

    return(ppi_logistic_ci(X, Y, Yhat, X_unlabeled, Yhat_unlabeled,

      alpha = alpha, optimizer_options = optimizer_options,

      alternative = alternative, lhat = lhat, coord = coord, w = w,

      w_unlabeled = w_unlabeled))
  }

  var_unlabeled <- cov(lhat * stats$grads_hat_unlabeled)

  cat("var_unlabeled:", var_unlabeled, "\n")

  var <- cov(stats$grads - lhat * stats$grads_hat)

  cat("var:", var, "\n")

  Sigma_hat <- stats$inv_hessian %*% (n / N * var_unlabeled + var) %*% stats$inv_hessian

  cat("Sigma_hat:", Sigma_hat, "\n")

  return(zconfint_generic(ppi_pointest, sqrt(diag(Sigma_hat) / n),

    alpha = alpha, alternative = alternative))
}

#=== END =======================================================================
