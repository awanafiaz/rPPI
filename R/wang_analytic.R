#===============================================================================
#
#  PROGRAM: wang_analytic.R
#
#  AUTHORS: Stephen Salerno (ssalerno@fredhutch.org)
#           Kentaro Hoffman (khoffm3@uw.edu)
#
#
#  PURPOSE: Implementation of various algorithms from Wang et al. (2020)
#
#           Methods for correcting inference based on outcomes predicted by
#           machine learning
#
#  INPUTS:  rel_form, inf_form, dat
#
#  OUTPUS:  Post-prediction inference functions for various target estimands:
#
#           1. wang_analy: Analytic Correction
#           2. wang_boot:  Bootstrap Correction
#           3. wang_mi: Multiple Imputation
#
#  Notes:   1. How do we want to handle the data argument? One argument with a
#              stacked dataset and a column that gives the set (tr, te, val)?
#              Or two dataset arguments, an unlabeled and labeled?
#
#  Updated: 2023-12-04
#
#===============================================================================

#=== ANALYTIC CORRECTION =======================================================

#' PPI Linear Regression using Wang et al. (2020) Analytic Correction
#'
#' @description
#' A short description...
#'
#' @details
#' Additional details...
#'
#' @param rel_form A formula defining the relationship model. This should be
#' of the form Y ~ Yhat, where Y is the name of the column corresponding to
#' the observed outcome in the labeled data and Yhat is the name of the column
#' corresponding to the predicted outcome in the labeled data.
#'
#' @param inf_form A formula defining the inference model. This should be of
#' the form Yhat ~ X, where Yhat is the name of the column corresponding to the
#' predicted outcome in the unlabeled data, and X generally corresponds to the
#' features of interest (e.g., X1 + X2).
#'
#' @param dat data in the form of the simdat function
#'
#' @returns description
#'
#' @examples
#'
#' rel_form <- Y ~ Yhat
#'
#' inf_form <- Yhat ~ X1
#'
#' dat <- simdat()
#'
#' wang_analytic(rel_form, inf_form, dat = dat)
#'
#' @export
#'
#' @import stats
#'

#-- WANG - ANALYTIC

wang_analytic <- function(rel_form, inf_form, dat) {

  #- 1. Estimate Prediction Model (Done in Data Step)

  #- 2. Estimate Relationship Model

  fit_rel <- lm(rel_form, data = dat[dat$set == "tst",])

  #- 3. Estimate Inference Model

  fit_inf <- lm(inf_form, data = dat[dat$set == "val",])

  #- 4. Coefficient Estimator

  X_val <- model.matrix(inf_form, data = dat[dat$set == "val",])

  est <- solve(crossprod(X_val)) %*% t(X_val) %*%

    (coef(fit_rel)[1] + coef(fit_rel)[2]*X_val %*% coef(fit_inf))

  #- 5. SE of Coefficient Estimator

  se <- sqrt(diag(solve(crossprod(X_val))*(sigma(fit_rel)^2 +

                                             (coef(fit_rel)[2]^2)*sigma(fit_inf)^2)))

  #- Output

  return(list(est = as.vector(est), se = as.vector(se)))
}

