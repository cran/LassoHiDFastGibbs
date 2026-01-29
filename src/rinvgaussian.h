 
#ifndef RINVGAUSSIAN_H
#define RINVGAUSSIAN_H
 
#include <RcppArmadillo.h>
 
arma::vec rrinvgauss(arma::vec vmu, arma::vec vlambda);
arma::vec rinvgauss_safe(arma::vec vmu, arma::vec vlambda);
arma::vec rinvgauss_slice(arma::vec vx, arma::vec vmu, arma::vec vlambda); 
  
#endif
