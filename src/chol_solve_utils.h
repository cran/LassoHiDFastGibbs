#ifndef FASTGIBBSSAMPLERS_CHOL_SOLVE_UTILS_H
#define FASTGIBBSSAMPLERS_CHOL_SOLVE_UTILS_H

#include <RcppArmadillo.h>

// IMPORTANT: mark these inline because they live in a header.
// Otherwise, every .cpp that includes this header will define them,
// causing "multiple definition" at link time.

inline arma::vec chol_solve_upper_vec(const arma::mat& mR, const arma::vec& vb)
{
  arma::vec vy = arma::solve(arma::trimatl(mR.t()), vb);
  arma::vec vx = arma::solve(arma::trimatu(mR), vy);
  return vx;
}

inline arma::mat chol_solve_upper_mat(const arma::mat& mR, const arma::mat& mB)
{
  arma::mat mY = arma::solve(arma::trimatl(mR.t()), mB);
  arma::mat mX = arma::solve(arma::trimatu(mR), mY);
  return mX;
}

inline arma::vec chol_solve_lower_vec(const arma::mat& mL, const arma::vec& vb)
{
  arma::vec vy = arma::solve(arma::trimatl(mL), vb);
  arma::vec vx = arma::solve(arma::trimatu(mL.t()), vy);
  return vx;
}

inline arma::mat chol_solve_lower_mat(const arma::mat& mL, const arma::mat& mB)
{
  arma::mat mY = arma::solve(arma::trimatl(mL), mB);
  arma::mat mX = arma::solve(arma::trimatu(mL.t()), mY);
  return mX;
}

#endif  // FASTGIBBSSAMPLERS_CHOL_SOLVE_UTILS_H
