#' Data generation for examples
#'
#' @param n vector of size 3 indicating the sample size in the training, labelled, and unlabelled data sets
#' @param beta1 first regression coefficient (or, regression coefficient of variable of interest for inference)
#  @inheritParams gam::gam
#'
#' @return A data frame containing 4 regressors, labelled outcome, predicted outcome and a character variable indicating which dat set the observartion belongs to (training, test, validation).
#'
#' @export
#'
#' @import stats gam
#'
#'
#' @examples
#' simdat(c(100, 100, 100), 1)
#'


simdat <- function(n = c(300, 300, 300), beta1 = 1) {

  X1 <- rnorm(sum(n), 1)
  X2 <- rnorm(sum(n), 1)
  X3 <- rnorm(sum(n), 1)
  X4 <- rnorm(sum(n), 2)

  Y <- c(beta1*X1 + 0.5*X2 + 3*smooth(X3) + 4*smooth(X4) + rnorm(sum(n)))

  set <- rep(c("trn", "tst", "val"), n)

  dat <- data.frame(X1, X2, X3, X4, Y, Yhat = NA, set)

  fit_gam <- gam::gam(Y ~ s(X1) + s(X2) + s(X3) + s(X4), data = dat[set == "trn",])

  dat[set == "tst", "Yhat"] <- predict(fit_gam, newdat = dat[set == "tst",])

  dat[set == "val", "Yhat"] <- predict(fit_gam, newdat = dat[set == "val",])

  return(dat)
}
