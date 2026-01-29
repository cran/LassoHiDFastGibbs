
/**
 * This function can be used for
 *  - Both n>p and p>n cases
 *  - Lasso and horseshoe penalties
 *  - 2 Block Gibbs with 3/1 blocks
 *  - nested Gibbs sampling
 */

#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include "rinvgaussian.h"
#include "slice_samplers.h"

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

// [[Rcpp::export]]
Rcpp::List lm_penalized_nested_gibbs(const arma::vec& vy,
                                     const arma::mat& mX,
                                     const Rcpp::String penalty_type,
                                     const double& a,
                                     const double& b,
                                     const double& u,
                                     const double& v,
                                     const int& nsamples,
                                     const double& lambda_init,
                                     const arma::vec& va_init,
                                     const int& verbose,
                                     const double& lower,
                                     const double& upper,
                                     const int s_beta=1,
                                     const int s_siglam=1)
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

  const int maxiter = nsamples;

  // Constant values
  const double a_til = a + 0.5*n;
  // const double u_til = u + 0.5*p;

  // Initialise storage of samples for MCMC
  arma::mat mBeta(maxiter*s_beta,p);
  arma::vec vsigma2(maxiter*s_beta*s_siglam);
  arma::vec vlambda2(maxiter*s_beta*s_siglam);

  // Initialize parameter values
  arma::vec vbeta;
  arma::vec va = va_init;
  double sigma2;
  double sigma;
  double lambda = lambda_init;
  double lambda2 = lambda*lambda;

  // Declare full conditional variables
  arma::vec vmu_til;
  arma::vec vsigma2_til;
  double b_til;
  // double v_til;

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

  int count_siglam = 0;
  int count_beta = 0;

  for (int i = 0; i < maxiter; ++i)
  {
    va_inv = 1/va; // Cost O(p)
    va_invsqrt = sqrt(va_inv); // Cost O(p)

    if (n >= p)
    {
      mZ = XTX.each_row() % va_invsqrt.as_row(); // Cost O(p^2)
      mZ = mZ.each_col() % va_invsqrt; // Cost O(p^2)
      mZ = symmatu(mZ); // Cost O(p^2)
      eig_sym(vs2, mV, mZ, "dc"); // Cost O(p^3)
      vs2.clamp(1.0E-8,1.0E32);
      vs = sqrt(vs2); // Cost O(p)
      vc = mV.t() * (va_invsqrt % XTy ); // Cost O(p^2)
    }
    else
    {
      XAinv = mX.each_row() % va_inv.t(); // Costs O(np)
      mZ = XAinv*mX.t(); // Costs O(pn^2) - To do run in parallel?
      mZ = symmatu(mZ);  // Costs O(n^2)
      try {
        eig_sym(vs2, mU, mZ, "dc"); // Costs O(n^3)

        // I think we need to do this otherwise the slice sampler breaks,
        // or at least breaks more easily.
        vs2.clamp(1.0E-8,1.0E32); // Costs O(n)

      } catch (const std::exception &e) {
        Rcpp::stop("Error in eig_sym: %s", e.what());
      }
      vs = sqrt(vs2); // Costs O(n)
      vc = vs % (mU.t()*vy); // Cost O(n^2)
    }

    vc2 = square(vc);   // Cost O(m)

    ////////////////////////////////////////////////////////////////////////////

    for (int j = 0; j < s_beta; ++j)
    {
      for (int k = 0; k < s_siglam; ++k)
      {
        double u_star = u + 0.5*min_np - 1.0;
        lambda2 = slice_orthogonal(lambda2, vs2, vc2, a_til, b, u_star, v, yTy, lower, upper);
        vlambda2[count_siglam] = lambda2;

        ////////////////////////////////////////////////////////////////////////////

        // Sampling from sigma2|data,lambda2,va ~ IG()
        vw = 1.0/(vs2 + lambda2); // Costs O(m)
        double sigma2_hat = (yTy - dot(vc2, vw))/n; // Costs O(m)
        b_til = b + 0.5*n*sigma2_hat;
        sigma2 = 1.0/randg(distr_param(a_til,1.0/b_til)); // Costs O(1)
        sigma = sqrt(sigma2);
        vsigma2[count_siglam] = sigma2;

        count_siglam++;
      }

      // Sampling from beta|rest~ N()
      if (n >= p) {
        vec vw_sqrt = sqrt(vw);
        vbeta = va_invsqrt % (mV * (vw % vc + sigma * vw_sqrt % randn(p))); // O(p^2)
      } else {
        //vmu_til = va_inv % (mX.t() * (mU * (vw % (mU.t() * vy) )));
        //vu = sqrt(sigma2/lambda2)*(va_invsqrt % randn(p)); // O(p)
        //vv = mX*vu + sigma*randn(n); // O(np)
        //vbeta = vmu_til + vu - va_inv % (mX.t() * (mU * (vw % (mU.t() * vv) ))); // O(np)

        vu = sqrt(sigma2/lambda2)*(va_invsqrt % randn(p)); // O(p)
        vv = vy - mX*vu + sqrt(sigma2)*randn(n);          // O(np)
        vbeta = vu + va_inv % (mX.t() * (mU * (vw % (mU.t() * vv) ))); // O(np)
      }
      mBeta.row(count_beta) = vbeta.as_row();
      count_beta++;

      ////////////////////////////////////////////////////////////////////////////

      // Sample from va|beta,sig^2,lam^2 ~ IG()

      if (pen_type==1) {
        // Lasso penalty
        vnu = sqrt(sigma2/lambda2)/abs(vbeta); // O(p)
        va = rinvgauss_safe(vnu, vone_p); // O(p)
      }

      if (pen_type==2) {
        // Horseshoe penalty
        arma::vec vc_vals = 0.5*lambda2*square(vbeta)/sigma2;
        for (arma::uword j = 0; j < p; ++j)
        {
          double b_val = sqrt(va[j]);
          b_val = sample_eta(b_val, vc_vals[j]);
          va[j] = b_val*b_val;
        }
        va.clamp(1.0E-12,1.0E12);
      }

      ////////////////////////////////////////////////////////////////////////////
    }

    if (verbose!=0) {
      if ((i%verbose)==0) {
        Rcout << "iter: " << i << " lambda2: " << lambda2 << " sigma2: " << sigma2 << "\n";
      }
    }

  }

  return List::create(_["mBeta"] = mBeta,
                      _["vsigma2"] = vsigma2,
                      _["vlambda2"] = vlambda2
  );
}



