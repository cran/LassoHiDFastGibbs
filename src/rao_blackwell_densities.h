
#ifndef RAO_BLACKWELL_DENSITIES_H
#define RAO_BLACKWELL_DENSITIES_H

#include <RcppArmadillo.h>

////////////////////////////////////////////////////////////////////////////////

arma::vec rb_normal(arma::vec vx, arma::vec vmu, arma::vec vsigma2);
arma::vec rb_invgamma(arma::vec vx, arma::vec va, arma::vec vb);
arma::vec rb_gamma(arma::vec vx, arma::vec va, arma::vec vb);
arma::vec rb_lasso(arma::vec vx, arma::vec va, arma::vec vb, arma::vec vc);

////////////////////////////////////////////////////////////////////////////////

#endif
