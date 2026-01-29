
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include "rinvgaussian.h"
#include "slice_samplers.h"
#include "chol_solve_utils.h"


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

using namespace Rcpp;
using namespace arma;
using namespace Eigen;
using namespace Numer;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

// arma::vec chol_solve_upper_vec(const arma::mat& mR, const arma::vec& vb)
// {
  // Step 1: Solve R^T * y = b (forward substitution)
//  arma::vec vy = arma::solve(arma::trimatl(mR.t()), vb);

  // Step 2: Solve R * x = y (backward substitution)
//  arma::vec vx = arma::solve(arma::trimatu(mR), vy);

//  return vx;
// }


// arma::mat chol_solve_upper_mat(const arma::mat& mR, const arma::mat& mB)
// {
  // Step 1: Solve R^T * y = b (forward substitution)
//  arma::mat mY = arma::solve(arma::trimatl(mR.t()), mB);

  // Step 2: Solve R * x = y (backward substitution)
//  arma::mat mX = arma::solve(arma::trimatu(mR), mY);

//  return mX;
// }


// arma::vec chol_solve_lower_vec(const arma::mat& mL, const arma::vec& vb)
// {
  // Forward substitution to solve L * y = b
//  arma::vec vy = arma::solve(arma::trimatl(mL), vb);

  // Backward substitution to solve L^T * x = y
//  arma::vec vx = arma::solve(arma::trimatu(mL.t()), vy);

//   return vx;
// }


// arma::mat chol_solve_lower_mat(const arma::mat& mL, const arma::mat& mB)
// {
  // Forward substitution to solve L * y = b
//   arma::mat mY = arma::solve(arma::trimatl(mL), mB);

  // Backward substitution to solve L^T * x = y
//   arma::mat mX = arma::solve(arma::trimatu(mL.t()), mY);

//   return mX;
// }


////////////////////////////////////////////////////////////////////////////////

/**
 * lm_penalized_4bg
 *
 * Sampling from auxiliary variables has been modified form Park and Casella
 * paper to facilitate comparison with other Gibbs samplers. This function also
 * samples from lambda2. It is a 4-block Gibbs sampler written for comparing to the two block algorithm
 *
 * Code is suitable for n > p (no modification is made to speed up code for p>n case).
 *
 * Inputs:
 *  vy - an n-vector of responses.
 *  mX - an n by p design matrix. Code assumes n>p
 *  lambda - a positive scalar
 *  sigma2_init - a starting value for sigma2. Shouldn't matter overly much what value is provided.
 *  a - prior shape parameter for inverse gamma prior for sigma2
 *  b - prior scale parameter for inverse gamma prior for sigma2
 *  u - prior shape parameter for inverse gamma prior for lambda2
 *  v - prior scale parameter for inverse gamma prior for lambda2
 *  nsamples - number of samples to run algorithm for.
 *  verbose - whether or not to  display sampling information
 *
 *  Outputs
 *  mBeta - an nsamples by p matrix of samples for vbeta
 *  vsigma2 - a nsamples vector of samples for sigma2
 *  mM - if rao_blackwellization an nsamples by p matrix of the mean parameters of vbeta|rest
 *  mV - if rao_blackwellization an nsamples by p matrix of the variance parameters of vbeta|rest
 *  a_til - The shape parameter sigma2|rest ~ IG - is a constant
 *  vb_til - if rao_blackwellization an nsamples vector of scale parameters of sigma2|rest ~ IG
 *  u_til - The shape parameter lambda2|rest ~ IG - is a constant
 *  vv_til - if rao_blackwellization an nsamples vector of scale parameters of lambda2|rest ~ IG
 */


// [[Rcpp::export]]
Rcpp::List lm_penalized_4bg(const arma::vec& vy,
                      const arma::mat& mX,
                      const Rcpp::String penalty_type,
                      const double& lambda_init,
                      const double& sigma2_init,
                      const double& a,
                      const double& b,
                      const double& u,
                      const double& v,
                      const int& nsamples,
                      const int& verbose=10000)
{
  //////////////////////////////////////////////////////////////////////////////

  int pen_type = 0;
  if (penalty_type=="lasso") {
    pen_type = 1;
  }
  if (penalty_type=="horseshoe") {
    pen_type = 2;
  }
  if (pen_type==0) {
    stop("Error: penalty_type not supported. Current options are lasso or horseshoe.");
  }

  //////////////////////////////////////////////////////////////////////////////

  // Initialize dimensions of problem
  arma::uword n = mX.n_rows;
  arma::uword p = mX.n_cols;

  // Summary statistics that can be calculated once at the beginning
  const arma::vec vone_n = arma::ones(n);
  const arma::vec vone_p = arma::ones(p);
  const arma::vec XTy = mX.t() * vy;
  const double yTy = dot(vy, vy);

  arma::mat XTX;
  arma::mat mI_n;

  if (n > p) {
    XTX =  mX.t() * mX;
  } else {
    mI_n = eye(n, n);
  }

  int maxiter = nsamples;

  // Initialise storage of samples
  arma::mat mBeta(maxiter,p);
  arma::vec vsigma2(maxiter);
  arma::vec vlambda2(maxiter);

  // Initialization
  arma::vec va = arma::ones(p);
  double sigma2 = sigma2_init;          // Set current value of sigma2
  double lambda2 = lambda_init*lambda_init;       // lambda2 is lambda squared - is constant
  const double a_til = a + 0.5*(n + p); // Shape of sigma2|rest ~ IG - is constant
  double b_til;                         // Scale of sigma2|rest ~ IG
  const double u_til = u + 0.5*p;       // Shape of lambda2|rest ~ IG - is constant
  double v_til;                         // Shape of lambda2|rest ~ IG - is constant

  arma::mat vmu_til;    // Mean of beta|rest ~ N
  arma::mat mSigma_til; // Covariance of beta|rest ~ N
  arma::vec vbeta;      // Current value of vbeta
  arma::mat mQ;         // An intermediary value
  arma::mat mR;
  arma::mat mR_inv;

  arma::mat mZ;
  arma::mat mU, mV;
  arma::vec vs;
  arma::vec vs2;
  arma::vec vu;
  arma::vec vv;
  arma::vec vc;
  arma::vec vc2;
  arma::vec vw;
  arma::vec vnu;

  arma::vec va_inv;
  arma::vec va_invsqrt;
  arma::mat XAinv;
  // bool okay = true;

  arma::vec vbeta2;
  arma::vec vsigma2_til;

  // Parameters of vd|rest ~ Inv-Gaussian
  const arma::vec vlambda = arma::ones(p); // Constant
  double Q;
  bool chol_okay;

  // Main loop
  for (int i = 0; i < maxiter; ++i)
  {
    // Sample from vbeta|rest
    if (n>=p) {
      chol_okay = chol(mR, XTX + diagmat(lambda2*va));
      if (chol_okay)
      {
        mR_inv = mR.i();
        vmu_til = mR_inv * (mR_inv.t() * XTy);

        // Sample from vbeta|sigma2,rest
        vbeta = vmu_til + sqrt(sigma2) * (mR_inv * randn(p));
      }
      else
      {
        stop("Cholesky factorisation of XTX + diagmat(lambda2*vd) failed!");
        // Note: Could attempt to recover more gracefully here
        // For example, attempt QR or SVD instead of Cholesky factorization.
      }

    } else {

      va_inv = 1/va;
      va_invsqrt = sqrt(va_inv);

      XAinv = mX.each_row() % va_inv.t();
      mQ = XAinv*mX.t() + lambda2*mI_n;

      // Sometimes a warning is given that mQ is not symmetric
      mQ = symmatu(mQ);

      // This approach is faster but more prone to numerical instability
      chol_okay = chol(mR, mQ);
      if (!chol_okay) {
        stop("Error: cholesky factorization failed! Try setting use_chol to false.");
      }

      vmu_til = XAinv.t() * chol_solve_upper_vec(mR, vy);

      // Sample from vbeta|rest
      vu = sqrt(sigma2/lambda2) * va_invsqrt % randn(p);
      vv = mX*vu + sqrt(sigma2)*randn(n);
      vbeta = vmu_til + vu -   (XAinv.t() * chol_solve_upper_vec(mR, vv));
      vbeta2 = vbeta % vbeta;
    }

    Q = dot(va, square(vbeta));

    ////////////////////////////////////////////////////////////////////////////

    // Sample from sigma2|rest
    double RSS;
    if (n>=p) {
      RSS =  yTy + dot(vbeta, XTX * vbeta - 2.0*XTy); // O(p^2)
    } else {
      RSS = sum(square(vy - mX*vbeta)); // O(np)
    }
    b_til = b + 0.5*(RSS + lambda2*Q);  // O(p)
    sigma2 = 1/randg(distr_param(a_til,1/b_til)); // O(p)

    ////////////////////////////////////////////////////////////////////////////

    // Sample from lambda2|rest
    v_til = v + 0.5*Q/sigma2; // O(p)
    lambda2 = randg(distr_param(u_til,1/v_til)); // O(p)

    ////////////////////////////////////////////////////////////////////////////

    // Sample from vd|rest
    if (pen_type==1) {
      // Lasso penalty
      vnu = sqrt(sigma2/lambda2)/abs(vbeta); // O(p)
      va = rinvgauss_safe(vnu, vone_p); // O(p)
    }

    if (pen_type==2) {
      // Horseshoe penalty
      vec vc_vals = 0.5*lambda2*square(vbeta)/sigma2;
      for (arma::uword j = 0; j < p; ++j)
      {
        double b_val = sqrt(va[j]);
        b_val = sample_eta(b_val, vc_vals[j]);
        va[j] = b_val*b_val;
      }
      va.clamp(1.0E-12,1.0E12);
    }

    ////////////////////////////////////////////////////////////////////////////

    if (verbose!=0) {
      if ((i%verbose)==0) {
        Rcout << "iter: " << i << " lambda2: " << lambda2 << " sigma2: " << sigma2 << "\n";
      }
    }

    // Storing samples
    mBeta.row(i) = vbeta.as_row();
    vsigma2[i] = sigma2;
    vlambda2[i] = lambda2;
  }

  return List::create(_["mBeta"] = mBeta,
                      _["vsigma2"] = vsigma2,
                      _["vlambda2"] = vlambda2);
}


////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
Rcpp::List lm_penalized_3bg(const arma::vec& vy,
                            const arma::mat& mX,
                            const Rcpp::String penalty_type,
                            const double& lambda,
                            const double& sigma2_init,
                            const double& a,
                            const double& b,
                            const double& u,
                            const double& v,
                            const int& nsamples,
                            const int& verbose=10000)
{
  //////////////////////////////////////////////////////////////////////////////

  int pen_type = 0;
  if (penalty_type=="lasso") {
    pen_type = 1;
  }
  if (penalty_type=="horseshoe") {
    pen_type = 2;
  }
  if (pen_type==0) {
    stop("Error: penalty_type not supported. Current options are lasso or horseshoe.");
  }

  //////////////////////////////////////////////////////////////////////////////

  // Initialize dimensions of problem
  arma::uword n = mX.n_rows;
  arma::uword p = mX.n_cols;

  // Summary statistics that can be calculated once at the beginning
  const arma::vec vone_n = arma::ones(n);
  const arma::vec vone_p = arma::ones(p);
  const arma::vec XTy = mX.t() * vy;
  const double yTy = dot(vy, vy);

  arma::mat XTX;
  arma::mat mI_n;

  if (n >= p) {
    XTX =  mX.t() * mX;
  } else {
    mI_n = eye(n, n);
  }

  int maxiter = nsamples;

  // Initialise storage of samples
  arma::mat mBeta(maxiter,p);
  arma::vec vsigma2(maxiter);
  arma::vec vlambda2(maxiter);

  // Initialization
  arma::vec va = arma::ones(p);
  double sigma2 = sigma2_init;          // Set current value of sigma2
  // double sigma = sqrt(sigma2);
  double lambda2 = lambda*lambda;       // lambda2 is lambda squared - is constant
  const double a_til = a + 0.5*n;       // Shape of sigma2|rest ~ IG - is constant
  double b_til;                         // Scale of sigma2|rest ~ IG
  const double u_til = u + 0.5*p;       // Shape of lambda2|rest ~ IG - is constant
  double v_til;                         // Shape of lambda2|rest ~ IG - is constant

  arma::mat vmu_til;    // Mean of beta|rest ~ N
  arma::mat mSigma_til; // Covariance of beta|rest ~ N
  arma::vec vbeta;      // Current value of vbeta
  arma::mat mQ;         // An intermediary value
  arma::mat mR;
  arma::mat mR_inv;

  arma::mat mZ;
  arma::mat mU, mV;
  arma::vec vs;
  arma::vec vs2;
  arma::vec vu;
  arma::vec vv;
  arma::vec vc;
  arma::vec vc2;
  arma::vec vw;
  arma::vec vnu;
  arma::mat XAinv;
  arma::vec va_inv;
  arma::vec va_invsqrt;

  arma::vec vbeta2;
  arma::vec vsigma2_til;

  // Parameters of va|rest ~ Inv-Gaussian
  const arma::vec vlambda = arma::ones(p); // Constant
  // double Q;
  bool chol_okay;

  // Main loop
  for (int i = 0; i < maxiter; ++i)
  {

    if (n>=p) {
      chol_okay = chol(mR, XTX + diagmat(lambda2*va));
      if (chol_okay)
      {
        mR_inv = mR.i();
        vmu_til = mR_inv * (mR_inv.t() * XTy);

        // Sample from sigma2|rest
        double RSS  =  yTy - dot(XTy, vmu_til); // O(p)
        b_til = b + 0.5*RSS;  // O(p)
        sigma2 = 1/randg(distr_param(a_til,1/b_til)); // O(p)
        // double sigma = sqrt(sigma2);

        // Sample from vbeta|sigma2,rest
        vbeta = vmu_til + sqrt(sigma2) * (mR_inv * randn(p));
      }
      else
      {
        stop("Cholesky factorisation of XTX + diagmat(lambda2*vd) failed!");
        // Note: Could attempt to recover more gracefully here
        // For example, attempt QR or SVD instead of Cholesky factorization.
      }

    } else {
      va_inv = 1/va;
      va_invsqrt = sqrt(va_inv);
      XAinv = mX.each_row() % va_inv.t();
      mQ = XAinv*mX.t() + lambda2*mI_n;

      // Sometimes a warning is given that mQ is not symmetric
      mQ = symmatu(mQ);

      // This approach is faster but more prone to numerical instability
      chol_okay = chol(mR, mQ);
      if (!chol_okay) {
        stop("Error: cholesky factorization failed! Try setting use_chol to false.");
      }

      vmu_til = XAinv.t() * chol_solve_upper_vec(mR, vy);

      // Sample from sigma2|rest
      double sigma2_hat = (yTy - dot(XTy, vmu_til))/n;
      b_til = b + 0.5*n*sigma2_hat;
      sigma2 = 1/randg(distr_param(a_til,1/b_til));

      // Sample from vbeta|rest
      vu = sqrt(sigma2/lambda2) * va_invsqrt % randn(p);
      vv = mX*vu + sqrt(sigma2)*randn(n);
      vbeta = vmu_til + vu -   (XAinv.t() * chol_solve_upper_vec(mR, vv));
      vbeta2 = vbeta % vbeta;
    }

    ////////////////////////////////////////////////////////////////////////////

    // Sample from lambda2|rest
    v_til = v + 0.5*dot(va, square(vbeta))/sigma2; // O(p)
    lambda2 = randg(distr_param(u_til,1/v_til)); // O(p)

    ////////////////////////////////////////////////////////////////////////////

    // Sample from vd|rest
    if (pen_type==1) {
      // Lasso penalty
      vnu = sqrt(sigma2/lambda2)/abs(vbeta); // O(p)
      va = rinvgauss_safe(vnu, vone_p); // O(p)
    }

    if (pen_type==2) {
      // Horseshoe penalty
      vec vc_vals = 0.5*lambda2*square(vbeta)/sigma2;
      for (arma::uword j = 0; j < p; ++j)
      {
        double b_val = sqrt(va[j]);
        b_val = sample_eta(b_val, vc_vals[j]);
        va[j] = b_val*b_val;
      }
      va.clamp(1.0E-12,1.0E12);
    }

    ////////////////////////////////////////////////////////////////////////////

    if (verbose!=0) {
      if ((i%verbose)==0) {
        Rcout << "iter: " << i << " lambda2: " << lambda2 << " sigma2: " << sigma2 << "\n";
      }
    }

    // Storing samples
    mBeta.row(i) = vbeta.as_row();
    vsigma2[i] = sigma2;
    vlambda2[i] = lambda2;
  }

  return List::create(_["mBeta"] = mBeta,
                      _["vsigma2"] = vsigma2,
                      _["vlambda2"] = vlambda2);
}

////////////////////////////////////////////////////////////////////////////////

