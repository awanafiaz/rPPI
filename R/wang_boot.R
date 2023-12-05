#===============================================================================
#
#  PROGRAM: wang_boot.R
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
#  INPUTS:  rel_form, inf_form, dat, n
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

#=== BOOTSTRAP CORRECTION =======================================================

#' PPI Linear Regression using Wang et al. (2020) Bootstrap Correction
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
#' @param nboot number of bootstraps
#'
#' @param rel_func relationship function model form (lm, gam, glm, etc)
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
#' nboot <- 100
#'
#' wang_boot(rel_form, inf_form, dat = dat, nboot)
#'
#' @export
#'
#' @import stats
#' @importFrom ranger ranger
#' @importFrom gam gam


#-- WANG - BOOTSTRAP

wang_boot <- function(rel_form, inf_form, dat, nboot = 100, rel_func = "lm") {

  #-- 1. Estimate Prediction Model (Done in Data Step)

  #-- 2. Estimate Relationship Model

  if (rel_func == "lm") {

    fit_rel <- lm(rel_form, data = dat[dat$set == "tst",])

  } else if (rel_func == "rf") {

    fit_rel <- ranger(rel_form, data = dat[dat$set == "tst",], keep.inbag = T)

  } else if (rel_func == "gam") {

    fit_rel <- gam(rel_form, data = dat[dat$set == "tst",])

  } else {

    stop("Currently only 'lm', 'rf', and 'gam' are supported")
  }

  #-- 3. Bootstrap

  set.seed(12345)

  dat_val <- dat[dat$set == "val",]

  n_val <- nrow(dat_val)

  inf_form_b <- reformulate(all.vars(inf_form)[-1], response = "Y_tilde_b")

  ests_b <- sapply(1:nboot, function(b) {

    #-   i. Sample Predicted Values and Covariates with Replacement

    idx_b <- sample(1:n_val, n_val, replace = T)

    dat_val_b <- dat_val[idx_b, ]

    #-  ii. Simulate Values from Relationship Model

    if (rel_func == "lm") {

      dat_val_b$Y_tilde_b <- rnorm(

        n_val, predict(fit_rel, dat_val_b), sigma(fit_rel))

    } else if (rel_func == "rf") {

      rel_preds <- predict(fit_rel, data = dat_val_b, type = "se")

      dat_val_b$Y_tilde_b <- rnorm(n_val, rel_preds$predictions, rel_preds$se)

    } else if (rel_func == "gam") {

      dat_val_b$Y_tilde_b <- rnorm(

        n_val, predict(fit_rel, dat_val_b), sigma(fit_rel))

    } else {

      stop("Currently only 'lm', 'rf', and 'gam' are supported")
    }

    #- iii. Fit Inference Model on Simulated Outcomes

    fit_inf_b <- lm(inf_form_b, data = dat_val_b)

    #-  iv. Extract Coefficient Estimator

    #-   v. Extract SE of Estimator

    return(summary(fit_inf_b)$coefficients[, 1:2])
  })

  #-- 4. Estimate Inference Model Coefficient

  est <- apply(ests_b[1:(nrow(ests_b)/2),], 1, median)

  #-- 5. Estimate Inference Model SE

  #- a. Parametric Bootstrap

  se_p <- apply(ests_b[(nrow(ests_b)/2 + 1):nrow(ests_b),], 1, median)

  #- b. Nonparametric Bootstrap

  se_n <- apply(ests_b[1:(nrow(ests_b)/2),], 1, sd)

  #-- Output

  return(list(est = est, se_par = se_p, se_npar = se_n))
}
