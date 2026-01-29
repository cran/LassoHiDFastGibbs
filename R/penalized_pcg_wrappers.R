#' Penalized PCG sampler: lambda2 collapsed over sigma2
#'
#' Partially-collapsed Gibbs (PCG) sampler for a Gaussian linear model with a
#' shrinkage prior/penalty on regression coefficients. This variant samples
#' \eqn{\lambda^2} using a collapsed step over \eqn{\sigma^2} (see implementation).
#'
#' @param vy Numeric response vector of length n.
#' @param mX Numeric design matrix of dimension n x p.
#' @param penalty_type Character string: \code{"lasso"} or \code{"horseshoe"}.
#' @param a,b Hyperparameters for the inverse-gamma prior on \eqn{\sigma^2}.
#' @param u,v Hyperparameters for the prior on \eqn{\lambda^2}.
#' @param nsamples Number of MCMC iterations.
#' @param lambda_init Initial value for \eqn{\lambda}.
#' @param sigma2_init Initial value for \eqn{\sigma^2}.
#' @param verbose Print progress every \code{verbose} iterations (0 = silent).
#'
#' @return A list with components:
#' \describe{
#'   \item{mBeta}{Matrix of beta draws (nsamples x p).}
#'   \item{vsigma2}{Vector of sigma^2 draws (length nsamples).}
#'   \item{vlambda2}{Vector of lambda^2 draws (length nsamples).}
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 40; p <- 6
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- c(1.2, 2, -1, 0.5, 0.75, 2.5)
#' y <- X %*% beta + rnorm(n)
#' out <- penalized_pcg_lambda2_sigma2(
#'   vy = y, mX = X, penalty_type = "lasso",
#'   a = 1, b = 1, u = 1, v = 1,
#'   nsamples = 200, lambda_init = 1, sigma2_init = 1,
#'   verbose = 0
#' )
#' str(out)
#'
#' @export
penalized_pcg_lambda2_sigma2 <- function(vy, mX, penalty_type = c("lasso", "horseshoe"),
                                         a, b, u, v, nsamples,
                                         lambda_init = 1, sigma2_init = 1,
                                         verbose = max(1L, floor(as.integer(nsamples) / 10))) {
  penalty_type <- match.arg(penalty_type)
  lm_penalized_pcg_lambda2_col_sigma2(
    vy = as.numeric(vy),
    mX = as.matrix(mX),
    penalty_type = penalty_type,
    a = a, b = b, u = u, v = v,
    nsamples = as.integer(nsamples),
    lambda_init = lambda_init,
    sigma2_init = sigma2_init,
    verbose = as.integer(verbose)
  )
}

#' Penalized PCG sampler: sigma2 collapsed over lambda2
#'
#' Partially-collapsed Gibbs (PCG) sampler variant that samples \eqn{\sigma^2} using a collapsed step over
#' \eqn{\lambda^2} (see implementation). Requires initial values for local scales
#' \code{va_init}; if omitted, it is set to a vector of ones.
#'
#' @inheritParams penalized_pcg_lambda2_sigma2
#' @param va_init Optional initial local-scale vector (length p). If \code{NULL},
#'   defaults to \code{rep(1, ncol(mX))}.
#' @param lower,upper Bounds used by the slice sampler.
#'
#' @return A list with components:
#' \describe{
#'   \item{mBeta}{Matrix of beta draws (nsamples x p).}
#'   \item{vsigma2}{Vector of sigma^2 draws (length nsamples).}
#'   \item{vlambda2}{Vector of lambda^2 draws (length nsamples).}
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 40; p <- 6
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- c(1.2, 2, -1, 0.5, 0.75, 2.5)
#' y <- X %*% beta + rnorm(n)
#'
#' out <- penalized_pcg_sigma2_lambda2(
#'   vy = y, mX = X, penalty_type = "lasso",
#'   a = 1, b = 1, u = 1, v = 1,
#'   nsamples = 200,
#'   lambda_init = 1, sigma2_init = 1,
#'   va_init = rep(1, p),
#'   verbose = 0
#' )
#'
#' str(out)

#'
#' @export
penalized_pcg_sigma2_lambda2 <- function(vy, mX, penalty_type = c("lasso", "horseshoe"),
                                         a, b, u, v, nsamples,
                                         lambda_init = 1, sigma2_init = 1,
                                         va_init = NULL,
                                         verbose = max(1L, floor(as.integer(nsamples) / 10)),
                                         lower = 1e-12, upper = 5000) {
  penalty_type <- match.arg(penalty_type)
  mX <- as.matrix(mX)
  if (is.null(va_init)) va_init <- rep(1, ncol(mX))

  lm_penalized_pcg_sigma2_col_lambda2(
    vy = as.numeric(vy),
    mX = mX,
    penalty_type = penalty_type,
    a = a, b = b, u = u, v = v,
    nsamples = as.integer(nsamples),
    lambda_init = lambda_init,
    va_init = as.numeric(va_init),
    verbose = as.integer(verbose),
    lower = lower,
    upper = upper
  )
}

#' Penalized PCG sampler: beta block, lambda2 collapsed over sigma2
#'
#' Partially-collapsed Gibbs (PCG) sampler variant that updates \eqn{\beta} in a dedicated block and samples
#' \eqn{\lambda^2} using a collapsed step over \eqn{\sigma^2}.
#'
#' @inheritParams penalized_pcg_lambda2_sigma2
#'
#' @return A list with components:
#' \describe{
#'   \item{mBeta}{Matrix of beta draws (nsamples x p).}
#'   \item{vsigma2}{Vector of sigma^2 draws (length nsamples).}
#'   \item{vlambda2}{Vector of lambda^2 draws (length nsamples).}
#' }
#'
#'
#' @examples
#' set.seed(1)
#' n <- 40; p <- 6
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- c(1.2, 2, -1, 0.5, 0.75, 2.5)
#' y <- X %*% beta + rnorm(n)
#'
#' out <- penalized_pcg_beta_sigma2(
#'   vy = y, mX = X, penalty_type = "horseshoe",
#'   a = 1, b = 1, u = 1, v = 1,
#'   nsamples = 200,
#'   lambda_init = 1, sigma2_init = 1,
#'   verbose = 0
#' )
#'
#' summary(out$mBeta)
#'
#' @export
penalized_pcg_beta_sigma2 <- function(vy, mX, penalty_type = c("lasso", "horseshoe"),
                                      a, b, u, v, nsamples,
                                      lambda_init = 1, sigma2_init = 1,
                                      verbose = max(1L, floor(as.integer(nsamples) / 10))) {
  penalty_type <- match.arg(penalty_type)
  lm_penalized_pcg_vbeta_col_sigma2(
    vy = as.numeric(vy),
    mX = as.matrix(mX),
    penalty_type = penalty_type,
    a = a, b = b, u = u, v = v,
    nsamples = as.integer(nsamples),
    lambda_init = lambda_init,
    sigma2_init = sigma2_init,
    verbose = as.integer(verbose)
  )
}

#' Penalized PCG sampler: sigma2 collapsed over beta
#'
#' Partially-collapsed Gibbs (PCG) sampler variant that samples \eqn{\sigma^2} using a collapsed step over
#' \eqn{\beta}. Requires initial values for local scales \code{va_init}.
#'
#' @inheritParams penalized_pcg_sigma2_lambda2
#'
#' @return A list with components:
#' \describe{
#'   \item{mBeta}{Matrix of beta draws (nsamples x p).}
#'   \item{vsigma2}{Vector of sigma^2 draws (length nsamples).}
#'   \item{vlambda2}{Vector of lambda^2 draws (length nsamples).}
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 40; p <- 6
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- c(1.2, 2, -1, 0.5, 0.75, 2.5)
#' y <- X %*% beta + rnorm(n)
#'
#' out <- penalized_pcg_sigma2_beta(
#'   vy = y, mX = X, penalty_type = "horseshoe",
#'   a = 1, b = 1, u = 1, v = 1,
#'   nsamples = 200,
#'   lambda_init = 1, sigma2_init = 1,
#'   va_init = rep(1, p),
#'   verbose = 0
#' )
#'
#' summary(out$vsigma2)
#'
#' @export
penalized_pcg_sigma2_beta <- function(vy, mX, penalty_type = c("lasso", "horseshoe"),
                                      a, b, u, v, nsamples,
                                      lambda_init = 1, sigma2_init = 1,
                                      va_init = NULL,
                                      verbose = max(1L, floor(as.integer(nsamples) / 10)),
                                      lower = 1e-12, upper = 5000) {
  penalty_type <- match.arg(penalty_type)
  mX <- as.matrix(mX)
  if (is.null(va_init)) va_init <- rep(1, ncol(mX))

  lm_penalized_pcg_sigma2_col_vbeta(
    vy = as.numeric(vy),
    mX = mX,
    penalty_type = penalty_type,
    a = a, b = b, u = u, v = v,
    nsamples = as.integer(nsamples),
    lambda_init = lambda_init,
    va_init = as.numeric(va_init),
    verbose = as.integer(verbose),
    lower = lower,
    upper = upper
  )
}

#' Bayesian lasso PCG sampler: lambda2 collapsed over local scales
#'
#' Lasso-specific Partially-collapsed Gibbs (PCG) variant with the local scales (va) collapsed in the
#' \eqn{\lambda^2} update.
#'
#' @inheritParams penalized_pcg_lambda2_sigma2
#'
#' @return A list with components:
#' \describe{
#'   \item{mBeta}{Matrix of beta draws (nsamples x p).}
#'   \item{vsigma2}{Vector of sigma^2 draws (length nsamples).}
#'   \item{vlambda2}{Vector of lambda^2 draws (length nsamples).}
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 40; p <- 6
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- c(1.2, 2, -1, 0.5, 0.75, 2.5)
#' y <- X %*% beta + rnorm(n)
#'
#' out <- blasso_pcg_lambda2_va(
#'   vy = y, mX = X,
#'   a = 1, b = 1, u = 1, v = 1,
#'   nsamples = 200,
#'   lambda_init = 1, sigma2_init = 1,
#'   verbose = 0
#' )
#'
#' summary(out$vlambda2)
#'
#' @export
blasso_pcg_lambda2_va <- function(vy, mX, a, b, u, v, nsamples,
                                  lambda_init = 1, sigma2_init = 1,
                                  verbose = max(1L, floor(as.integer(nsamples) / 10))) {
  blasso_pcg_lambda2_col_va(
    vy = as.numeric(vy),
    mX = as.matrix(mX),
    a = a, b = b, u = u, v = v,
    nsamples = as.integer(nsamples),
    lambda_init = lambda_init,
    sigma2_init = sigma2_init,
    verbose = as.integer(verbose)
  )
}

#' Bayesian lasso PCG sampler: sigma2 collapsed over local scales
#'
#' Lasso-specific Partially-collapsed Gibbs (PCG) variant with the local scales (va) collapsed in the
#' \eqn{\sigma^2} update.
#'
#' @inheritParams penalized_pcg_sigma2_lambda2
#'
#' @return A list with components:
#' \describe{
#'   \item{mBeta}{Matrix of beta draws (nsamples x p).}
#'   \item{vsigma2}{Vector of sigma^2 draws (length nsamples).}
#'   \item{vlambda2}{Vector of lambda^2 draws (length nsamples).}
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 40; p <- 6
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- c(1.2, 2, -1, 0.5, 0.75, 2.5)
#' y <- X %*% beta + rnorm(n)
#'
#' out <- blasso_pcg_sigma2_va(
#'   vy = y, mX = X,
#'   a = 1, b = 1, u = 1, v = 1,
#'   nsamples = 200,
#'   lambda_init = 1, sigma2_init = 1,
#'   va_init = rep(1, p),
#'   verbose = 0
#' )
#'
#' summary(out$vsigma2)
#'
#' @export
blasso_pcg_sigma2_va <- function(vy, mX, a, b, u, v, nsamples,
                                 lambda_init = 1, sigma2_init = 1,
                                 va_init = NULL,
                                 verbose = max(1L, floor(as.integer(nsamples) / 10)),
                                 lower = 1e-12, upper = 5000) {
  mX <- as.matrix(mX)
  if (is.null(va_init)) va_init <- rep(1, ncol(mX))

  lm_penalized_pcg_sigma2_col_va(
    vy = as.numeric(vy),
    mX = mX,
    a = a, b = b, u = u, v = v,
    nsamples = as.integer(nsamples),
    lambda_init = lambda_init,
    sigma2_init = sigma2_init,
    va_init = as.numeric(va_init),
    verbose = as.integer(verbose),
    lower = lower,
    upper = upper
  )
}
