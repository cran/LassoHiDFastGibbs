#' Penalized nested Gibbs sampler for Bayesian linear regression
#'
#' Runs the nested Gibbs sampler for a Gaussian linear model
#' \eqn{y = X\beta + \epsilon} with either a lasso or horseshoe penalty
#' (shrinkage prior) on \eqn{\beta}. The algorithm supports both \eqn{n \ge p}
#' and \eqn{p > n} regimes.
#'
#' @param vy Numeric response vector of length \eqn{n}.
#' @param mX Numeric design matrix of dimension \eqn{n \times p}.
#' @param penalty_type Character string: either \code{"lasso"} or \code{"horseshoe"}.
#' @param a,b Hyperparameters for the inverse-gamma prior on \eqn{\sigma^2}.
#' @param u,v Hyperparameters for the prior on \eqn{\lambda^2}.
#' @param nsamples Integer number of outer MCMC iterations.
#' @param lambda_init Initial value for \eqn{\lambda}.
#' @param va_init Optional initial values for the local shrinkage parameters
#'   (length \eqn{p}). If \code{NULL}, a vector of ones is used.
#' @param verbose Print progress every \code{verbose} iterations (0 = silent).
#' @param lower,upper Bounds for the slice sampler used for \eqn{\lambda^2}.
#' @param s_beta Integer: number of inner updates of \eqn{\beta} per outer iteration.
#' @param s_siglam Integer: number of inner updates of \eqn{(\sigma^2,\lambda^2)}
#'   per \eqn{\beta} update.
#'
#' @return A list with components:
#' \describe{
#'   \item{mBeta}{Matrix of sampled \eqn{\beta} draws (rows correspond to stored draws).}
#'   \item{vsigma2}{Vector of sampled \eqn{\sigma^2} draws.}
#'   \item{vlambda2}{Vector of sampled \eqn{\lambda^2} draws.}
#' }
#'
#' \code{lm_penalized_nested_gibbs()}.
#'
#' @examples
#' set.seed(1)
#' n <- 50; p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' y <- rnorm(n)
#'
#' out <- penalized_nested_Gibbs(
#'   vy = y, mX = X,
#'   penalty_type = "lasso",
#'   a = 1, b = 1, u = 1, v = 1,
#'   nsamples = 200,
#'   lambda_init = 1,
#'   va_init = NULL,
#'   verbose = 0,
#'   lower = 1e-12,
#'   upper = 5000,
#'   s_beta = 1,
#'   s_siglam = 1
#' )
#' str(out)
#'
#' @export
penalized_nested_Gibbs <- function(
    vy,
    mX,
    penalty_type = c("lasso", "horseshoe"),
    a, b, u, v,
    nsamples,
    lambda_init = 1,
    va_init = NULL,
    verbose = max(1L, floor(as.integer(nsamples) / 10)),
    lower = 1e-12,
    upper = 5000,
    s_beta = 1L,
    s_siglam = 1L
) {
  penalty_type <- match.arg(penalty_type)
  mX <- as.matrix(mX)
  if (is.null(va_init)) va_init <- rep(1, ncol(mX))

  lm_penalized_nested_gibbs(
    vy = as.numeric(vy),
    mX = mX,
    penalty_type = penalty_type,
    a = a, b = b, u = u, v = v,
    nsamples = as.integer(nsamples),
    lambda_init = lambda_init,
    va_init = as.numeric(va_init),
    verbose = as.integer(verbose),
    lower = lower,
    upper = upper,
    s_beta = as.integer(s_beta),
    s_siglam = as.integer(s_siglam)
  )
}
