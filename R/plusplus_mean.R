#=== MEAN ESTIMATION ===========================================================

#' ...Need title...
#'
#' @description
#' Computes the prediction-powered point estimate of the p-dimensional mean.
#'
#' @details
#' `[ADZ23] <https://arxiv.org/abs/2311.01453>`__ A. N. Angelopoulos, J. C.
#' Duchi, and T. Zrnic. PPI++: Efficient Prediction Powered Inference.
#' arxiv:2311.01453, 2023.
#'
#' @param Y (ndarray): Gold-standard labels.
#'
#' @param Yhat (ndarray): Predictions corresponding to the gold-standard
#' labels.
#'
#' @param Yhat_unlabeled (ndarray): Predictions corresponding to the unlabeled
#' data.
#'
#' @param lhat (float, optional): Power-tuning parameter (see
#' `[ADZ23] <https://arxiv.org/abs/2311.01453>`__). The default value `None`
#' will estimate the optimal value from data. Setting `lhat=1` recovers PPI
#' with no power tuning, and setting `lhat=0` recovers the classical point
#' estimate.
#'
#' @param coord (int, optional): Coordinate for which to optimize `lhat`. If
#' `None`, it optimizes the total variance over all coordinates. Must be in
#' {1, ..., d} where d is the dimension of the estimand.
#'
#' @param w (ndarray, optional): Sample weights for the labeled data set.
#' Defaults to all ones vector.
#'
#' @param w_unlabeled (ndarray, optional): Sample weights for the unlabeled
#' data set. Defaults to all ones vector.
#'
#' @returns float or ndarray: Prediction-powered point estimate of the mean.
#'
#' @examples
#'
#' #need examples
#'
#' @import stats
#'
#' @export

ppi_mean_pointestimate <- function(Y, Yhat, Yhat_unlabeled,

  lhat = NULL, coord = NULL, w = NULL, w_unlabeled = NULL) {

  #- Compute Dimensions of Inputs

  n <- ifelse(is.null(dim(Y)), length(Y), nrow(Y))
  N <- ifelse(is.null(dim(Yhat_unlabeled)), length(Yhat_unlabeled), nrow(Yhat_unlabeled))
  d <- if (length(dim(Yhat)) > 1) dim(Yhat)[2] else 1

  #- Set Default Weights if Not Provided

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

  # If lhat is not provided, estimate it

  if (is.null(lhat)) {

    ppi_pointest <- mean(w_unlabeled * Yhat_unlabeled) + mean(w * (Y - Yhat))

    grads <- w * (Y - ppi_pointest)

    grads_hat <- w * (Yhat - ppi_pointest)

    grads_hat_unlabeled <- w_unlabeled * (Yhat_unlabeled - ppi_pointest)

    inv_hessian <- diag(d)

    lhat <- calc_lhat_glm(grads, grads_hat, grads_hat_unlabeled, inv_hessian, coord, clip = T)

    return(ppi_mean_pointestimate(Y, Yhat, Yhat_unlabeled, lhat = lhat, coord = coord, w = w, w_unlabeled = w_unlabeled))

  } else {

    return(mean(w_unlabeled * lhat * Yhat_unlabeled, na.rm = T) + mean(w * (Y - lhat * Yhat), na.rm = T))
  }
}

#' ...Need title...
#'
#' @description
#' Computes the prediction-powered confidence interval for a d-dimensional mean.
#'
#' @details
#' `[ADZ23] <https://arxiv.org/abs/2311.01453>`__ A. N. Angelopoulos, J. C.
#' Duchi, and T. Zrnic. PPI++: Efficient Prediction Powered Inference.
#' arxiv:2311.01453, 2023.
#'
#' @param Y (ndarray): Gold-standard labels.
#'
#' @param Yhat (ndarray): Predictions corresponding to the gold-standard
#' labels.
#'
#' @param Yhat_unlabeled (ndarray): Predictions corresponding to the unlabeled
#' data.
#'
#' @param alpha (float, optional): Error level; the confidence interval will target a coverage of 1 - alpha. Must be in (0, 1).
#'
#' @param alternative (str, optional): Alternative hypothesis, either 'two-sided', 'larger' or 'smaller'.
#'
#' @param lhat (float, optional): Power-tuning parameter (see
#' `[ADZ23] <https://arxiv.org/abs/2311.01453>`__). The default value `None`
#' will estimate the optimal value from data. Setting `lhat=1` recovers PPI
#' with no power tuning, and setting `lhat=0` recovers the classical point
#' estimate.
#'
#' @param coord (int, optional): Coordinate for which to optimize `lhat`. If
#' `None`, it optimizes the total variance over all coordinates. Must be in
#' {1, ..., d} where d is the dimension of the estimand.
#'
#' @param w (ndarray, optional): Sample weights for the labeled data set.
#' Defaults to all ones vector.
#'
#' @param w_unlabeled (ndarray, optional): Sample weights for the unlabeled
#' data set. Defaults to all ones vector.
#'
#' @returns tuple: Lower and upper bounds of the prediction-powered confidence
#' interval for the mean.
#'
#' @examples
#'
#' #need examples
#'
#' @import stats
#'
#' @export

ppi_mean_ci <- function(Y, Yhat, Yhat_unlabeled,

  alpha = 0.05, alternative = "two-sided", lhat = NULL, coord = NULL, w = NULL,

  w_unlabeled = NULL) {

  #- Compute Dimensions of Inputs

  n <- ifelse(is.null(dim(Y)), length(Y), nrow(Y))
  N <- ifelse(is.null(dim(Yhat_unlabeled)), length(Yhat_unlabeled), nrow(Yhat_unlabeled))
  d <- if (length(dim(Yhat)) > 1) dim(Yhat)[2] else 1

  #- Set Default Weights if Not Provided

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

  # If lhat is not provided, estimate it

  if (is.null(lhat)) {

    ppi_pointest <- ppi_mean_pointestimate(Y, Yhat, Yhat_unlabeled, lhat = 1, w = w, w_unlabeled = w_unlabeled)

    grads <- w * (Y - ppi_pointest)

    grads_hat <- w * (Yhat - ppi_pointest)

    grads_hat_unlabeled <- w_unlabeled * (Yhat_unlabeled - ppi_pointest)

    inv_hessian <- diag(d)

    lhat <- calc_lhat_glm(grads, grads_hat, grads_hat_unlabeled, inv_hessian, coord = NULL, clip = T)

    return(ppi_mean_ci(Y, Yhat, Yhat_unlabeled, lhat = lhat, coord = coord, w = w, w_unlabeled = w_unlabeled))
  }

  ppi_pointest <- ppi_mean_pointestimate(Y, Yhat, Yhat_unlabeled, lhat = lhat, coord = coord, w = w, w_unlabeled = w_unlabeled)

  imputed_std <- sd(w_unlabeled * (lhat * Yhat_unlabeled)) * sqrt((N - 1) / N) / sqrt(N)

  rectifier_std <- sd(w * (Y - lhat * Yhat)) * sqrt((n - 1) / n) / sqrt(n)

  ### QC

  cat("ppi_pointest:\n", ppi_pointest, "\n")

  cat("imputed_std:\n", imputed_std, "\n")

  cat("rectifier_std:\n", rectifier_std, "\n")

  ###

  return(zconfint_generic(ppi_pointest, sqrt(imputed_std^2 + rectifier_std^2), alpha, alternative))
}

#' ...Need title...
#'
#' @description
#' Computes the prediction-powered p-value for a 1D mean.
#'
#' @details
#' `[ADZ23] <https://arxiv.org/abs/2311.01453>`__ A. N. Angelopoulos, J. C.
#' Duchi, and T. Zrnic. PPI++: Efficient Prediction Powered Inference.
#' arxiv:2311.01453, 2023.
#'
#' @param Y (ndarray): Gold-standard labels.
#'
#' @param Yhat (ndarray): Predictions corresponding to the gold-standard
#' labels.
#'
#' @param Yhat_unlabeled (ndarray): Predictions corresponding to the unlabeled
#' data.
#'
#' @param null (float): Value of the null hypothesis to be tested.
#'
#' @param alternative (str, optional): Alternative hypothesis, either 'two-sided', 'larger' or 'smaller'.
#'
#' @param lhat (float, optional): Power-tuning parameter (see
#' `[ADZ23] <https://arxiv.org/abs/2311.01453>`__). The default value `None`
#' will estimate the optimal value from data. Setting `lhat=1` recovers PPI
#' with no power tuning, and setting `lhat=0` recovers the classical point
#' estimate.
#'
#' @param coord (int, optional): Coordinate for which to optimize `lhat`. If
#' `None`, it optimizes the total variance over all coordinates. Must be in
#' {1, ..., d} where d is the dimension of the estimand.
#'
#' @param w (ndarray, optional): Sample weights for the labeled data set.
#' Defaults to all ones vector.
#'
#' @param w_unlabeled (ndarray, optional): Sample weights for the unlabeled
#' data set. Defaults to all ones vector.
#'
#' @returns float or ndarray: Prediction-powered p-value for the mean.
#'
#' @examples
#'
#' #need examples
#'
#' @import stats
#'
#' @export

ppi_mean_pval <- function(Y, Yhat, Yhat_unlabeled,

  null = 0, alternative = "two-sided", lhat = NULL, coord = NULL, w = NULL,

  w_unlabeled = NULL) {

  #- Compute Dimensions of Inputs

  n <- ifelse(is.null(dim(Y)), length(Y), nrow(Y))
  N <- ifelse(is.null(dim(Yhat_unlabeled)), length(Yhat_unlabeled), nrow(Yhat_unlabeled))

  #- Set Default Weights if Not Provided

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

  # If lhat is not provided, estimate it

  if (is.null(lhat)) {

    if (length(dim(Y)) > 1 && dim(Y)[2] > 1) {

      lhat <- 1

    } else {

      ppi_pointest <- mean(w_unlabeled * Yhat_unlabeled) + mean(w * (Y - Yhat))

      grads <- w * (Y - ppi_pointest)

      grads_hat <- w * (Yhat - ppi_pointest)

      grads_hat_unlabeled <- w_unlabeled * (Yhat_unlabeled - ppi_pointest)

      inv_hessian <- matrix(1, nrow=1, ncol=1)

      lhat <- calc_lhat_glm(grads, grads_hat, grads_hat_unlabeled, inv_hessian, coord = NULL)
    }
  }

  return(rectified_p_value(mean(w * Y - lhat * w * Yhat), sd(w * Y - lhat * w * Yhat) * sqrt((n - 1) / n) / sqrt(n),

    mean(w_unlabeled * lhat * Yhat_unlabeled), sd(w_unlabeled * lhat * Yhat_unlabeled) * sqrt((N - 1) / N) / sqrt(N), null, alternative))
}

#=== END =======================================================================
