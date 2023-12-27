#===============================================================================
#
#  PROGRAM: predpi_ols.R
#
#  AUTHORS: Kentaro Hoffman (khoffm3@uw.edu)
#           Stephen Salerno (ssalerno@fredhutch.org)
#
#           ^ UPDATE THESE
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
#  NOTES:   1. NEED TO ADD REFERENCES
#
#  Updated: 2023-12-27
#
#===============================================================================

#=== ORDINARY LEAST SQUARES ====================================================

#' PPI Linear Regression using Angelopoulos et al. (2023) Analytic Form
#'
#' @description
#' A short description...
#'
#' @details
#' Additional details...
#'
#' @param rec_form A formula defining the rectifier model. This should be of     # Redo these to conform to new conventions
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
#' angelopoulos_original(rec_form, inf_form, dat = dat)         ## Update with wrapper function or keep method-specific examples?
#'
#' @import stats
#'
#' @export

predpi_ols <- function(X_l, Y_l, f_l, X_u, f_u, n, p, N, ...) {

  #- 1. Prediction-Powered Estimator

  theta_tilde_f <- solve(crossprod(X_u)) %*% t(X_u) %*% f_u

  delta_hat_f   <- solve(crossprod(X_l)) %*% t(X_l) %*% (f_l - Y_l)

  theta_hat_pp  <- theta_tilde_f - delta_hat_f

  #- 2. Meat and Bread for Imputed Estimate

  Sigma_tilde <- crossprod(X_u) / N

  M_tilde <- sapply(1:N, function(i) {

    (c(f_u[i] - crossprod(X_u[i,], theta_tilde_f)))^2 *

      tcrossprod(X_u[i,])}) |>

    rowMeans() |> matrix(nrow = p)

  iSigma_tilde <- solve(Sigma_tilde)

  #- 3. Sandwich Variance Estimator for Imputed Estimate

  V_tilde <- iSigma_tilde %*% M_tilde %*% iSigma_tilde

  #- 4. Meat and Bread for Empirical Rectifier

  Sigma <- crossprod(X_l) / n

  M <- sapply(1:n, function(i) {

    (c(f_l[i] - Y_l[i] - crossprod(X_l[i,], delta_hat_f)))^2 *

      tcrossprod(X_l[i,])}) |>

    rowMeans() |> matrix(nrow = p)

  iSigma <- solve(Sigma)

  #- 5. Sandwich Variance Estimator for Empirical Rectifier

  V <- iSigma %*% M %*% iSigma

  #- 6. Standard Error Estimates

  se <- sqrt(diag(V) / n + diag(V_tilde) / N)

  #- Output

  return(list(est = as.vector(theta_hat_pp), se = as.vector(se)))
}

#=== END =======================================================================
