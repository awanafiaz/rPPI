#===============================================================================
#
#  PROGRAM: ppi.R
#
#  AUTHORS: Awan Afiaz (aafiaz@uw.edu)
#           Kentaro Hoffman (khoffm3@uw.edu)
#           Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: NEED
#
#  INPUTS:  N/A
#
#  OUTPUTS: NEED
#
#  Notes:   1. Propose using one formula argument, as the Wang methods are the
#              only ones to deviate from the "rectifier" framework, but the
#              relationship model formula could be created internally and we
#              could note this in the details.
#
#           2. Propose one data argument where the unlabeled/labeled data are
#              "stacked" -- need to decide if we want an explicit argument for
#              a column that indicates which rows belong to which set, or can
#              we assume that any rows with missing Y are unlabeled? OR have
#              two data arguments, where if the second is null, we assume the
#              first is the stacked data.
#
#           3. Propose having a "method" argument to signify which authors'
#              method to use and a "estimand" argument (maybe better name?) to
#              signify what model (e.g., lm, glm, ...). Suggest we concatenate
#              these two arguments to match to the appropriate helper function,
#              all of which will have the same naming convention:
#              method_estimand
#
#           4. All other arguments that relate to all the methods (e.g., alpha),
#              or all other method-specific arguments and will have defaults.
#
#           5. Propose all of the other methods take matrix arguments, to
#              minimize additional parsing after wrapper function assertions
#              and data pre-processing.
#
#  Updated: 2023-12-14
#
#===============================================================================

#=== MAIN PPI FUNCTION =========================================================

#' Valid and Efficient Post-Prediction Inference (PPI) using State-of-the-Art
#' Methods
#'
#' @description
#' A short description...
#'
#' @details
#' Additional details...
#'
#' @param formula description ... of the form Y - f ~ X1 + X2 + ...
#'
#' @param estimand description
#'
#' @param method description
#'
#' @param data description
#'
#' @param set ... column name for lab versus unlab...
#'
#' @param seed description
#'
#' @param alpha description
#'
#' @param alternative description
#'
#' @param ... Further arguments passed to or from specific methods. See details.
#'
#' @returns description
#'
#' @examples
#' # example code
#'
#'
#' @import stats
#'
#' @export

ppi <- function(formula, estimand, method, data, set, seed = NULL,

  alpha = 0.05, alternative = "two-sided", ...) {

  #--- CHECKS & ASSERTIONS -----------------------------------------------------

  #-- CHECK FOR DATA

  if (missing(data)) data <- environment(formula)

  #-- CHECK FOR VALID FAMILY & METHOD

  if (!(method %in% c("postpi", "leanpostpi", "predpi", "plusplus"))) {          # Check these

    stop(paste("'method' must be one of",

      "c('postpi', 'leanpostpi', 'predpi', 'plusplus').",                        # Check these

      "See the 'Details' section of the documentation for more information."))
  }

  if (!(estimand %in%

    c("mean", "ols", "binomial", "multinomial", "quantile"))) {

    stop(paste("'estimand' must be one of",

      "c('mean', 'ols', 'binomial', 'multinomial', 'quantile')",

      "See the 'Details' section of the documentation for more information."))
  }

  #-- SET SEED

  if (!is.null(seed)) set.seed(seed)

  #--- DATA --------------------------------------------------------------------

  #-- LABELED DATA

  X_l <- model.matrix(formula, data = data[data$set == "tst",])                  # Fix column argument

  Y_l <- data[data$set == "tst", all.vars(formula)[1]] |> matrix(ncol = 1)       # Fix column argument

  f_l <- data[data$set == "tst", all.vars(formula)[2]] |> matrix(ncol = 1)       # Fix column argument

  #-- UNLABELED DATA

  X_u <- model.matrix(formula, data = data[data$set == "val",])                  # Fix column argument

  f_u <- data[data$set == "val", all.vars(formula)[2]] |> matrix(ncol = 1)       # Fix column argument

  #-- DIMENSIONS

  n <- nrow(X_l)

  p <- ncol(X_l)

  N <- nrow(X_u)

  #--- METHOD ------------------------------------------------------------------

  func <- get(paste(method, estimand, sep = "_"))

  fit <- func(X_l, Y_l, f_l, X_u, f_u, n, p, N, ...)

  ci  <- zconfint_generic(fit$est, fit$se, alpha, alternative)

  #--- RETURN ------------------------------------------------------------------

  obj <- list(est = fit$est, se = fit$se, ci = ci, ...)

  class(obj) <- c("ppi", fit$class)

  return(obj)
}

#=== END =======================================================================
