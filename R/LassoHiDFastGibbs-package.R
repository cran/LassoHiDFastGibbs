#' @title Fast High-Dimensional Gibbs Samplers for Bayesian Lasso Regression
#'
#' @description Provides fast and scalable Gibbs sampling algorithms for
#' Bayesian Lasso regression model in high-dimensional
#' settings. The package implements efficient partially collapsed
#' and nested Gibbs samplers for Bayesian Lasso, with a focus on
#' computational efficiency when the number of predictors is large
#' relative to the sample size.
#' @name LassoHiDFastGibbs
#' @docType package
#' @useDynLib LassoHiDFastGibbs, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
