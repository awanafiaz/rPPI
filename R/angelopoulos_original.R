#===============================================================================
#
#  PROGRAM: angelopoulos.R (PPI original)
#
#  AUTHORS: Kentaro Hoffman (khoffm3@uw.edu)
#           Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: Implementation of various algorithms from Angelopoulos et al. (2023)
#
#           Prediction-Powered Inference
#
#  INPUTS:  N/A
#
#  OUTPUS:  Prediction-powered inference functions for various target estimands: # SS: Wrapper function around these?
#
#           1. angelopoulos: Linear Regression                                   # SS: Better naming conventions?
#
#  Notes:   1. How do we want to handle the data argument? One argument with a
#              stacked dataset and a column that gives the set (tr, te, val)?
#              Or two dataset arguments, an unlabeled and labeled?
#
#  Updated: 2023-12-02
#
#===============================================================================

#=== LINEAR REGRESSION =========================================================

#' PPI Linear Regression using Angelopoulos et al. (2023) Analytic Form
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
#' @returns A list of outputs: estimate of inference model parameters and corresponding standard error based on both parametric and non-parametric bootstrap methods.
#'
#' @examples
#'
#' rec_form <- Y - Yhat ~ X1
#'
#' inf_form <- Yhat ~ X1
#'
#' dat <- simdat()
#'
#' angelopoulos_original(rec_form, inf_form, dat = dat)
#'
#' @export
#'
#' @import stats
#'


#-- ANGELOPOULOS

# SS: taking out the defaults so this breaks with no input. Can force assertions later
# SS: See note, should we have two data arguments to avoid slicing?

angelopoulos_original <- function(rec_form, inf_form, dat) {

  X <- model.matrix(rec_form, data = dat[dat$set == "tst",])

  Y <- dat[dat$set == "tst", all.vars(rec_form)[1]]

  f <- dat[dat$set == "tst", all.vars(rec_form)[2]]

  X_tilde <- model.matrix(inf_form, data = dat[dat$set == "val",])

  f_tilde <- dat[dat$set == "val", all.vars(inf_form)[1]]

  #- 1. Prediction-Powered Estimator

  theta_tilde_f <- solve(crossprod(X_tilde))%*% t(X_tilde) %*% f_tilde

  delta_hat_f <- solve(crossprod(X)) %*% t(X) %*% (f - Y)

  theta_hat_pp <- theta_tilde_f - delta_hat_f

  #- 2. Meat and Bread for Imputed Estimate

  Sigma_tilde <- crossprod(X_tilde) / nrow(X_tilde)

  M_tilde <- sapply(1:nrow(X_tilde), function(i) {

    (c(f_tilde[i] - crossprod(X_tilde[i,], theta_tilde_f)))^2 *

      tcrossprod(X_tilde[i,])}) |>

    rowMeans() |> matrix(nrow = ncol(X_tilde))

  iSigma_tilde <- solve(Sigma_tilde)

  #- 3. Sandwich Variance Estimator for Imputed Estimate

  V_tilde <- iSigma_tilde %*% M_tilde %*% iSigma_tilde

  #- 4. Meat and Bread for Empirical Rectifier

  Sigma <- crossprod(X) / nrow(X)

  M <- sapply(1:nrow(X), function(i) {

    (c(f[i] - Y[i] - crossprod(X[i,], delta_hat_f)))^2 * tcrossprod(X[i,])}) |>

    rowMeans() |> matrix(nrow = ncol(X))

  iSigma <- solve(Sigma)

  #- 5. Sandwich Variance Estimator for Empirical Rectifier

  V <- iSigma %*% M %*% iSigma

  #- 6. Standard Error Estimate

  se <- sqrt((diag(V) / nrow(X)) + (diag(V_tilde) / nrow(X_tilde)))

  #- Output

  return(list(est = as.vector(theta_hat_pp), se = as.vector(se)))
}

#=== END =======================================================================
