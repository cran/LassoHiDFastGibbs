
#ifndef ZETA_H
#define ZETA_H

#include <RcppArmadillo.h>
 
double zetaOneLentz_c(double x, double tol, int maxiter);
double zetaOne_c(double x);
//arma::vec zeta_c(int k, arma::vec x);
double zeta_c(int k, double x);

#endif
