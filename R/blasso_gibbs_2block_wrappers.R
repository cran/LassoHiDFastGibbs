#' Bayesian lasso Gibbs sampler: 2-block (beta–sigma2) variant
#'
#' @description
#' Implements a two-block Gibbs sampler for the Bayesian lasso regression model
#' in which the regression coefficients are updated jointly with the noise variance
#'  \eqn{\sigma^2} in one block, while the global shrinkage parameter and
#' local shrinkage parameters are updated conditionally in separate steps.
#'
#' @param vy Numeric response vector of length n.
#' @param mX Numeric design matrix of dimension n x p.
#' @param a,b Hyperparameters for the inverse-gamma prior on sigma^2.
#' @param u,v Hyperparameters for the prior on lambda^2.
#' @param nsamples Integer number of MCMC iterations.
#' @param lambda_init Initial value for lambda.
#' @param sigma2_init Initial value for sigma^2.
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
#' n <- 30; p <- 6
#' X <- matrix(rnorm(n * p), n, p)
#' y <- rnorm(n)
#' out <- blasso_gibbs_2block_bs(
#'   vy = y, mX = X,
#'   a = 1, b = 1, u = 1, v = 1,
#'   nsamples = 200,
#'   lambda_init = 1, sigma2_init = 1,
#'   verbose = 0
#' )
#' str(out)
#'
#' @export
blasso_gibbs_2block_bs <- function(vy, mX,
                                   a, b, u, v,
                                   nsamples,
                                   lambda_init = 1,
                                   sigma2_init = 1,
                                   verbose = max(1L, floor(nsamples / 5))) {
  blasso_2bg_vbeta_sigma2(
    vy = vy, mX = mX,
    a = a, b = b, u = u, v = v,
    nsamples = as.integer(nsamples),
    lambda_init = lambda_init,
    sigma2_init = sigma2_init,
    verbose = as.integer(verbose)
  )
}

#' Bayesian lasso Gibbs sampler: 2-block (beta–lambda2) variant
#'
#' @description
#' Implements a two-block Gibbs sampler for the Bayesian lasso regression model
#' in which the regression coefficients are updated jointly with the global
#' shrinkage parameter \eqn{\lambda^2} in one block, while the noise variance and
#' local shrinkage parameters are updated conditionally in separate steps.
#'
#' @param vy Numeric response vector of length n.
#' @param mX Numeric design matrix of dimension n x p.
#' @param a,b Hyperparameters for the inverse-gamma prior on sigma^2.
#' @param u,v Hyperparameters for the prior on lambda^2.
#' @param nsamples Integer number of MCMC iterations.
#' @param lambda_init Initial value for lambda.
#' @param sigma2_init Initial value for sigma^2.
#' @param va_init Optional initial values for local shrinkage parameters (length p).
#' @param verbose Print progress every \code{verbose} iterations (0 = silent).
#' @param lower,upper Bounds used by the slice sampler for lambda^2.
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
#' n <- 30; p <- 6
#' X <- matrix(rnorm(n * p), n, p)
#' y <- rnorm(n)
#' out <- blasso_gibbs_2block_bl(
#'   vy = y, mX = X,
#'   a = 1, b = 1, u = 1, v = 1,
#'   nsamples = 200,
#'   lambda_init = 1, sigma2_init = 1,
#'   verbose = 0
#' )
#' str(out)
#'
#' @export
blasso_gibbs_2block_bl <- function(vy, mX,
                                   a, b, u, v,
                                   nsamples,
                                   lambda_init = 1,
                                   sigma2_init = 1,
                                   va_init = NULL,
                                   verbose = max(1L, floor(nsamples / 5)),
                                   lower = 1.0e-12,
                                   upper = 5000) {
  mX <- as.matrix(mX)
  p <- ncol(mX)

  if (is.null(va_init)) {
    va_init <- rep(1, p)
  }

  blasso_2bg_vbeta_lambda2(
    vy = vy, mX = mX,
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
