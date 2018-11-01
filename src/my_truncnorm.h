#ifndef _carfima_MY_TRUNCNORM_H
#define _carfima_MY_TRUNCNORM_H
#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Implementation of Truncated Normal Distribution Functions : from MSM package
// 1. for 'rtruncnorm', n=1 for ease of use.
// 2. for 'dtruncnorm', it also cares about a single number case.
double rtruncnorm(double mean, double sd, double lbd, double rbd);
double dtruncnorm(double x, double mean, double sd, double lower, double upper);


#endif
