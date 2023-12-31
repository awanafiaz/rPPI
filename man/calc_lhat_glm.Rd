% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{calc_lhat_glm}
\alias{calc_lhat_glm}
\title{Estimate PPI++ Power Tuning Parameter}
\usage{
calc_lhat_glm(
  grads,
  grads_hat,
  grads_hat_unlabeled,
  inv_hessian,
  coord = NULL,
  clip = FALSE
)
}
\arguments{
\item{grads}{(matrix): n x p matrix gradient of the loss function with
respect to the parameter evaluated at the labeled data.}

\item{grads_hat}{(matrix): n x p matrix gradient of the loss function with
respect to the model parameter evaluated using predictions on the labeled
data.}

\item{grads_hat_unlabeled}{(matrix): N x p matrix gradient of the loss
function with respect to the parameter evaluated using predictions on the
unlabeled data.}

\item{inv_hessian}{(matrix): p x p matrix inverse of the Hessian of the
loss function with respect to the parameter.}

\item{coord}{(int, optional): Coordinate for which to optimize \code{lhat}.
If \code{None}, it optimizes the total variance over all coordinates.
Must be in {1, ..., d} where d is the shape of the estimand.}

\item{clip}{(bool, optional): Whether to clip the value of lhat to be
non-negative. Defaults to \code{False}.}
}
\value{
(float): Optimal value of \code{lhat} in [0,1].
}
\description{
Calculates the optimal value of lhat for the prediction-powered confidence
interval for GLMs.
}
\examples{

dat <- simdat()

form <- Y - Yhat ~ X1

X_l <- model.matrix(form, data = dat[dat$set == "tst",])

Y_l <- dat[dat$set == "tst", all.vars(form)[1]] |> matrix(ncol = 1)

f_l <- dat[dat$set == "tst", all.vars(form)[2]] |> matrix(ncol = 1)

X_u <- model.matrix(form, data = dat[dat$set == "val",])

f_u <- dat[dat$set == "val", all.vars(form)[2]] |> matrix(ncol = 1)

n <- nrow(X_l)

p <- ncol(X_l)

N <- nrow(X_u)

est <- plusplus_ols_est(X_l, Y_l, f_l, X_u, f_u, n, p, N)

w_l <- rep(1, n)

w_u <- rep(1, N)

stats <- ols_get_stats(est, X_l, Y_l, f_l, X_u, f_u, w_l, w_u, use_u = TRUE)

calc_lhat_glm(stats$grads, stats$grads_hat, stats$grads_hat_unlabeled,
              stats$inv_hessian, coord = NULL, clip = FALSE)

}
