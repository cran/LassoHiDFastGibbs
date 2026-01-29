
/**
 * This code is to calculate Rao-Blackwellised estimates of various densities.
 */

#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include "lasso_distribution.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

using namespace Rcpp;
using namespace arma;
using namespace Eigen;
using namespace Numer;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

/**
 * Calculate the Rao-Blackwellization density estimate for the normal case
 *
 * Inputs
 * vx - set of points used to evaluate density
 * vmu - nsamples vector of means
 * vsigma2 - nsamples vector of variances
 *
 * Output
 * vf - Rao-Blackwellization density estimate
 */


arma::vec rb_normal(arma::vec vx, arma::vec vmu, arma::vec vsigma2)
{
  arma::uword N = vx.n_elem;
  arma::uword nsamples = vmu.n_elem;

  arma::vec vsigma = sqrt(vsigma2);

  arma::vec vf = arma::zeros(N);
  for (arma::uword i = 0; i < nsamples; ++i) {
    for (arma::uword j = 0; j < N; ++j) {
      vf(j) += normpdf(vx(j), vmu[i], vsigma[i]);
    }
  }
  vf = vf / nsamples;
  return vf;
}


////////////////////////////////////////////////////////////////////////////////

/**
 * Calculate the Rao-Blackwellization density estimate for the inverse gamma case
 *
 * Inputs
 * vx - set of points used to evaluate density
 * vmu - nsamples vector of means
 * vsigma2 - nsamples vector of variances
 *
 * Output
 * vf - Rao-Blackwellization density estimate
 */


arma::vec rb_invgamma(arma::vec vx, arma::vec va, arma::vec vb)
{
  arma::uword N = vx.n_elem;
  arma::uword nsamples = vb.n_elem;

  arma::vec vf = arma::zeros(N);
  for (arma::uword i = 0; i < nsamples; ++i) {
    for (arma::uword j = 0; j < N; ++j) {
      vf(j) += exp(va[i]*log(vb[i]) - lgamma(va[i]) - (va[i] + 1)*log(vx[j]) - vb[i]/vx[j]);
    }
  }
  vf = vf / nsamples;
  return vf;
}

////////////////////////////////////////////////////////////////////////////////

/**
 * Calculate the Rao-Blackwellization density estimate for the gamma case
 *
 * Inputs
 * vx - set of points used to evaluate density
 * vmu - nsamples vector of means
 * vsigma2 - nsamples vector of variances
 *
 * Output
 * vf - Rao-Blackwellization density estimate
 */


arma::vec rb_gamma(arma::vec vx, arma::vec va, arma::vec vb)
{
  arma::uword N = vx.n_elem;
  arma::uword nsamples = vb.n_elem;

  arma::vec vf = arma::zeros(N);
  for (arma::uword i = 0; i < nsamples; ++i) {
    for (arma::uword j = 0; j < N; ++j) {
      vf(j) += exp(va[i]*log(vb[i]) - lgamma(va[i]) + (va[i] - 1)*log(vx[j]) - vb[i]*vx[j]);
    }
  }
  vf = vf / nsamples;
  return vf;
}

////////////////////////////////////////////////////////////////////////////////


/**
 * Calculate the Rao-Blackwellization density estimate for the lasso distribution case
 *
 * Inputs
 * vx - set of points used to evaluate density
 * va - scale parameter
 * vb - location parameter
 * vc - penalty parameter
 *
 * Output
 * vf - Rao-Blackwellization density estimate
 */


arma::vec rb_lasso(arma::vec vx, arma::vec va, arma::vec vb, arma::vec vc)
{
  arma::uword N = vx.n_elem;
  arma::uword nsamples = vb.n_elem;

  arma::vec vf = arma::zeros(N);
  for (arma::uword i = 0; i < nsamples; ++i) {
    vf += dlasso_c_v2(vx, va[i], vb[i], vc[i], false);
  }
  vf = vf / nsamples;
  return vf;
}

////////////////////////////////////////////////////////////////////////////////
