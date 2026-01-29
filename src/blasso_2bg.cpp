
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include "rinvgaussian.h"
#include "slice_samplers.h"
#include "chol_solve_utils.h"


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

using namespace Rcpp;
using namespace arma;
using namespace Eigen;
using namespace Numer;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

// arma::vec chol_solve_upper_vec(const arma::mat& mR, const arma::vec& vb)
// {
  // Step 1: Solve R^T * y = b (forward substitution)
//   arma::vec vy = arma::solve(arma::trimatl(mR.t()), vb);

  // Step 2: Solve R * x = y (backward substitution)
//  arma::vec vx = arma::solve(arma::trimatu(mR), vy);

//  return vx;
// }


// arma::mat chol_solve_upper_mat(const arma::mat& mR, const arma::mat& mB)
// {
  // Step 1: Solve R^T * y = b (forward substitution)
//   arma::mat mY = arma::solve(arma::trimatl(mR.t()), mB);

  // Step 2: Solve R * x = y (backward substitution)
//  arma::mat mX = arma::solve(arma::trimatu(mR), mY);

//  return mX;
// }


// arma::vec chol_solve_lower_vec(const arma::mat& mL, const arma::vec& vb)
// {
  // Forward substitution to solve L * y = b
//   arma::vec vy = arma::solve(arma::trimatl(mL), vb);

  // Backward substitution to solve L^T * x = y
//  arma::vec vx = arma::solve(arma::trimatu(mL.t()), vy);

//  return vx;
// }


// arma::mat chol_solve_lower_mat(const arma::mat& mL, const arma::mat& mB)
// {
  // Forward substitution to solve L * y = b
//  arma::mat mY = arma::solve(arma::trimatl(mL), mB);

  // Backward substitution to solve L^T * x = y
//  arma::mat mX = arma::solve(arma::trimatu(mL.t()), mY);

//  return mX;
// }

////////////////////////////////////////////////////////////////////////////////


// [[Rcpp::export]]
Rcpp::List blasso_2bg_vbeta_sigma2(const arma::vec& vy,
                                const arma::mat& mX,
                                const double& a,
                                const double& b,
                                const double& u,
                                const double& v,
                                const int& nsamples,
                                const double& lambda_init,
                                const double& sigma2_init,
                                const int& verbose=1000)
{
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
  double sigma = sqrt(sigma2);
  double lambda2 = lambda_init*lambda_init;  // lambda2 is lambda squared - is constant
  const double a_til = a + 0.5*n;       // Shape of sigma2|rest ~ IG - is constant
  double b_til;                         // Scale of sigma2|rest ~ IG
  const double u_til = u + 0.5*p;       // Shape of lambda2|rest ~ IG - is constant
  // double v_til = v;                         // Shape of lambda2|rest ~ IG - is constant

  arma::mat vmu_til;    // Mean of beta|rest ~ N
  arma::mat mSigma_til; // Covariance of beta|rest ~ N
  arma::vec vbeta;      // Current value of vbeta
  arma::mat mQ;         // An intermediary value
  arma::mat mR;
  arma::mat mR_inv;

  arma::vec vbeta2;
  arma::vec vsigma2_til;

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

  bool chol_okay;

  // Parameters of vd|rest ~ Inv-Gaussian
  const arma::vec vlambda = arma::ones(p); // Constant
  // double Q;

  // Main loop
  for (int i = 0; i < maxiter; ++i)
  {
    ////////////////////////////////////////////////////////////////////////////

    // Block 1. Sample from vbeta,sigma2|rest
    va_inv = 1/va;
    va_invsqrt = sqrt(va_inv);

    if (n>=p) {
      bool okay = chol(mR, XTX + diagmat(lambda2*va));
      if (okay)
      {
        mR_inv = mR.i();
        vmu_til = mR_inv * (mR_inv.t() * XTy);

        // Sample from sigma2|rest
        double RSS  =  yTy - dot(XTy, vmu_til); // O(p)
        b_til = b + 0.5*RSS;  // O(p)
        sigma2 = 1/randg(distr_param(a_til,1/b_til)); // O(p)
        sigma = sqrt(sigma2);

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

      XAinv = mX.each_row() % va_inv.t();
      mQ = XAinv*mX.t() + lambda2*mI_n;

      // Sometimes a warning is given that mQ is not symmetric
      mQ = symmatu(mQ);

      // This approach is faster but more prone to numerical instability
      chol_okay = chol(mR, mQ);
      if (!chol_okay) {
        stop("Error: cholesky factorization failed!");
      }

      vmu_til = XAinv.t() * chol_solve_upper_vec(mR, vy);

      // Sample from sigma2|rest
      double sigma2_hat = (yTy - dot(XTy, vmu_til))/n;
      b_til = b + 0.5*n*sigma2_hat;
      sigma2 = 1/randg(distr_param(a_til,1/b_til));
      sigma = sqrt(sigma2);

      // Sample from vbeta|rest
      vu = sqrt(sigma2/lambda2) * va_invsqrt % randn(p);
      vv = mX*vu + sigma*randn(n);
      vbeta = vmu_til + vu - (XAinv.t() * chol_solve_upper_vec(mR, vv));
    }

    vbeta2 = square(vbeta);

    ////////////////////////////////////////////////////////////////////////////

    // Block 2. Sample from lambda2,va|rest

    // Sample from lambda2|rest
    double a_val = u_til - 1;
    double b_val = v;
    double c_val = sum(abs(vbeta))/sigma;

    lambda2 = slice_sampler_precision(lambda2, a_val, b_val, c_val);

    ////////////////////////////////////////////////////////////////////////////

    // Sample from va|lambda2,rest
    vnu = sqrt(sigma2/lambda2)/abs(vbeta);
    va = rinvgauss_safe(vnu, vlambda);

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
Rcpp::List blasso_2bg_vbeta_lambda2(const arma::vec& vy,
                                     const arma::mat& mX,
                                     const double& a,
                                     const double& b,
                                     const double& u,
                                     const double& v,
                                     const int& nsamples,
                                     const double& lambda_init,
                                     const double& sigma2_init,
                                     const arma::vec& va_init,
                                     const int& verbose,
                                     const double& lower,
                                     const double& upper)
{
  // Initialize dimensions of problem
  arma::uword n = mX.n_rows;
  arma::uword p = mX.n_cols;

  arma::vec np(2);
  np(0) = n;
  np(1) = p;
  arma::uword min_np = np.min();
  // arma::uword max_np = np.max();

  // Pre-calcuate algorithm constants
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

  // Constant values
  const double a_til = a + 0.5*n;
  // const double u_til = u + 0.5*p;

  const int maxiter = nsamples;

  // Initialise storage of samples for MCMC
  arma::mat mBeta(maxiter,p);
  arma::vec vsigma2(maxiter);
  arma::vec vlambda2(maxiter);


  // Initialize parameter values
  arma::vec vbeta;
  arma::vec va = va_init;
  double sigma2 = sigma2_init;
  double sigma = sqrt(sigma2);
  double lambda = lambda_init;
  double lambda2 = lambda*lambda;

  // Declare full conditional variables
  arma::vec vmu_til(p);
  arma::vec vsigma2_til(p);
  // double b_til = b;
  // double v_til = v;

  arma::mat mZ;
  arma::mat mU, mV;
  arma::vec vs;
  arma::vec vs2;
  arma::vec vc;
  arma::vec vc2;
  arma::vec vw;
  arma::vec vu;
  arma::vec vv;
  arma::vec vnu;

  arma::vec va_inv;
  arma::vec va_invsqrt;
  arma::mat XAinv;
  double sigma2_hat;

  // bool calc_svd = true;

  for (int i = 0; i < maxiter; ++i)
  {
    va_inv = 1/va; // Cost O(p)
    va_invsqrt = sqrt(va_inv); // Cost O(p)

    if (n>=p) {
      mZ = XTX.each_row() % va_invsqrt.as_row(); // Cost O(p^2)
      mZ = mZ.each_col() % va_invsqrt; // Cost O(p^2)
      mZ = symmatu(mZ); // Cost O(p^2)
      eig_sym(vs2, mV, mZ, "dc"); // Cost O(p^3)
      vs2.clamp(1.0E-8,1.0E32);
      vs = sqrt(vs2); // Cost O(p)
      vc = mV.t() * (va_invsqrt % XTy ); // Cost O(p^2)
    } else {

      XAinv = mX.each_row() % va_inv.t();
      mZ = XAinv*mX.t();
      mZ = symmatu(mZ);
      try {
        eig_sym(vs2, mU, mZ, "dc"); // Faster than svd_econ

        // I think we need to do this otherwise the slice sampler breaks,
        // or at least breaks more easily.
        vs2.clamp(1.0E-8,1.0E32);

      } catch (const std::exception &e) {
        Rcpp::stop("Error in eig_sym: %s", e.what());
      }
      vs = sqrt(vs2);
      vc = vs % (mU.t()*vy); // Cost O(n^2)
    }


    vc2 = square(vc);      // Cost O(m)
    double u_star = u + 0.5*min_np - 1.0;
    lambda2 = slice_orthogonal(lambda2, vs2, vc2, a_til, b, u_star, v, yTy, lower, upper);

    ////////////////////////////////////////////////////////////////////////////

    vw = 1.0/(vs2 + lambda2); // Costs O(m)
    sigma2_hat = (yTy - dot(vc2, vw))/n; // Costs O(m)

    ////////////////////////////////////////////////////////////////////////////

    // Sampling from vbeta|data,lambda2,sigma2,va ~ N()

    if (n>=p) {
      vbeta = va_invsqrt % (mV * (vw % vc + sqrt(sigma2)* sqrt(vw) % randn(p))); // O(p^2)
    } else {
      vu = sqrt(sigma2/lambda2)*(va_invsqrt % randn(p)); // O(p)
      vv = (vy - mX*vu) + sigma*randn(n);                // O(np)
      vbeta = vu + va_inv % (mX.t() * (mU * (vw % (mU.t() * vv) ))); // O(np)
    }

    ////////////////////////////////////////////////////////////////////////////

    // Slice sample from sigma2|rest
    lambda = sqrt(lambda2);
    double a_val = (a + 0.5*(n+p)-1);
    double b_val = b + 0.5*n*sigma2_hat;
    double c_val = lambda*sum(abs(vbeta));

    // Slice sampler for tau=1/sigma2 and then invert.
    double tau = 1/sigma2;
    tau = slice_sampler_precision(tau, a_val, b_val, c_val);
    sigma2 = 1/tau;
    sigma = sqrt(sigma2);

    ////////////////////////////////////////////////////////////////////////////

    // Sample from va|beta,sig^2,lam^2 ~ IG()
    vnu = sqrt(sigma2/lambda2)/abs(vbeta); // O(p)
    va = rinvgauss_safe(vnu, vone_p); // O(p)

    ////////////////////////////////////////////////////////////////////////////

    if (verbose!=0) {
      if ((i%verbose)==0) {
        Rcout << "iter: " << i << " lambda2: " << lambda2 << " sigma2: " << sigma2 << "\n";
      }
    }

    ////////////////////////////////////////////////////////////////////////////

    // Store MCMC samples
    mBeta.row(i) = vbeta.as_row();
    vsigma2[i] = sigma2;
    vlambda2[i] = lambda2;
  }

  //clock.stop("clock_times");

  return List::create(_["mBeta"] = mBeta,
                      _["vsigma2"] = vsigma2,
                      _["vlambda2"] = vlambda2);
}

