
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include "rinvgaussian.h"
#include "slice_samplers.h"
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

// [[Rcpp::export]]
Rcpp::List lasso_gibbs_modified_hans(const arma::vec& vy,
                                     const arma::mat& mX,
                                     const double& a,
                                     const double& b,
                                     const double& u,
                                     const double& v,
                                     const int& nsamples,
                                     const arma::vec& vbeta_init,
                                     const double& lambda_init,
                                     const double& sigma2_init,
                                     int verbose,
                                     bool tune_lambda2,
                                     bool rao_blackwellization)
{
  // Initialize dimensions of the problem
  arma::uword n = mX.n_rows;
  arma::uword p = mX.n_cols;

  const int maxiter = nsamples;

  // Initialize storage of samples for MCMC
  arma::mat mBeta(maxiter,p);
  arma::vec vsigma2(maxiter);
  arma::vec vlambda2(maxiter);

  // Declare variables for Rao-Blackwellization
  arma::mat mA;
  arma::mat mB;
  arma::mat mC;

  // Initialize variables for Rao-Blackwellization
  if (rao_blackwellization) {
    mA = zeros(maxiter,p);
    mB = zeros(maxiter,p);
    mC = zeros(maxiter,p);
  }

  // Initialize problem constants
  const arma::vec vone_n = arma::ones(n);
  const arma::vec vone_p = arma::ones(p);
  arma::vec dgXTX  = (mX%mX).t() * vone_n; // Diagonal values of XTX
  arma::mat XTX;
  if (n>p) {
    XTX = mX.t() * mX;
  }
  arma::vec XTy = mX.t() * vy;
  double yTy = dot(vy, vy);

  // Set the current values of the parameters
  arma::vec vd = arma::ones(p);
  arma::vec vnu = arma::ones(p);

  // Assign initial values
  arma::vec vbeta = vbeta_init;
  double sigma2 = sigma2_init;
  double sigma = sqrt(sigma2);
  double lambda = lambda_init;
  double lambda2 = lambda*lambda;

  // Constant values
  const double a_til = a + 0.5*(n + p);
  const double u_til = u + 0.5*p;

  for (int i = 0; i < maxiter; ++i)
  {
    arma::vec va_vals = dgXTX/sigma2;
    arma::vec vb_vals = zeros(p);
    arma::vec vc_vals = vone_p*sqrt(lambda2/sigma2);

    arma::vec XTy_hat;
    arma::vec vy_hat;
    arma::vec vy_hat_mj;
    arma::vec vx_j;

    arma::vec vu = randu(p);
    if (n > p)
    {
      XTy_hat = XTX * vbeta;
      for (arma::uword j=0; j < p; ++j) {
        vx_j = XTX.col(j);
        XTy_hat = XTy_hat - vx_j*vbeta[j]; // This might not be exactly right
        vb_vals[j] = (XTy[j] - XTy_hat[j])/sigma2;
        vbeta[j] =  qlasso_fast_c_v2(vu[j], va_vals[j], vb_vals[j], vc_vals[j]);
        XTy_hat = XTy_hat + vx_j*vbeta[j];
      }
    }
    else
    {
      vy_hat = mX*vbeta;
      for (arma::uword j=0; j < p; ++j) {
        vx_j = mX.col(j);
        vy_hat_mj = vy_hat - vx_j*vbeta[j];
        vb_vals[j] = (XTy[j] -  dot(vx_j, vy_hat_mj))/sigma2;
        vbeta[j] = qlasso_fast_c_v2(vu[j], va_vals[j], vb_vals[j], vc_vals[j]);
        vy_hat = vy_hat_mj + vx_j*vbeta[j];
      }
    }

    // Store sufficient statistics for Rao-Blackwellization
    if (rao_blackwellization) {
      mA.row(i) = va_vals.as_row();
      mB.row(i) = vb_vals.as_row();
      mC.row(i) = vc_vals.as_row();
    }

    ////////////////////////////////////////////////////////////////////////////

    // Slice from lambda2|rest
    double sum_abs_vbeta = sum(abs(vbeta));
    double RSS;
    if (n>p) {
      RSS = yTy - 2*dot(vbeta, XTy) + dot(vbeta, XTX*vbeta); // O(p^2)
    } else {
      RSS = sum(pow(vy - vy_hat,2.0)); // O(n)
    }

    // Slice sample from sigma2|rest
    lambda = sqrt(lambda2);
    double a_val = (a_til-1);
    double b_val = b + 0.5*RSS;
    double c_val = lambda*sum_abs_vbeta;

    // Slice sampler for tau=1/sigma2 and then invert.
    double tau = 1/sigma2;
    tau = slice_sampler_precision(tau, a_val, b_val, c_val);
    sigma2 = 1/tau;
    sigma = sqrt(sigma2);

    ////////////////////////////////////////////////////////////////////////////

    // Slice sampler for lambda2
    if (tune_lambda2)
    {
      a_val = u_til - 1;
      b_val = v;
      c_val = sum_abs_vbeta/sigma;

      lambda2 = slice_sampler_precision(lambda2, a_val, b_val, c_val);
      lambda  = sqrt(lambda2);
    }
    ////////////////////////////////////////////////////////////////////////////

    if (verbose!=0) {
      if ((i%verbose)==0) {
        Rcout << "iter: " << i << " lambda2: " << lambda2 << " sigma2: " << sigma2 << "\n";
      }
    }

    // Store MCMC samples
    mBeta.row(i) = vbeta.as_row();
    vsigma2[i] = sigma2;
    vlambda2[i] = lambda2;
  }

  return List::create(_["mBeta"] = mBeta,
                      _["vsigma2"] = vsigma2,
                      _["vlambda2"] = vlambda2,
                      _["mA"] = mA,
                      _["mB"] = mB,
                      _["mC"] = mC);
}



////////////////////////////////////////////////////////////////////////////////
