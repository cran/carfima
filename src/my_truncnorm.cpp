// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "my_truncnorm.h"

using namespace Rcpp;
using namespace arma;

/*
 * Implementation of Truncated Normal Distribution Functions : from MSM package
 * 1. rtruncnorm
 * 2. dtruncnorm
 * 
 * ptruncnorm and qtruncnorm are not implemented. Write them if you need those.
 * Note that in http://dirk.eddelbuettel.com/code/rcpp/html/Rmath_8h_source.html#l00035,
 * the element 'lt' is for 0 (lower.tail=FALSE) and 1 (lower.tail=TRUE)
 *             'lg' is for 0 (log.p = FALSE)    and 1 (log.p = TRUE)
 */

// 1. rtruncnorm
// 3-case argument
//    (1) Rejection Sampling with exponential proposal. use if lower >> mean
// [[Rcpp::export]]
double rtruncnorm(double mean, double sd, double lbd, double rbd){
  RNGScope rngScope;
  double output  = 0.0;
  double lower   = (lbd-mean)/sd;
  double upper   = (rbd-mean)/sd;
  bool rindicate = true;
  
  double a, z, u, rho;
  if ((lower>=0)&&(upper>lower+2.0*std::sqrt(std::exp(1.0))/(lower+std::sqrt(std::pow(lower,2)+4))*std::exp(2.0*lower - lower*std::sqrt(std::pow(lower,2.0)+4.0)/4))){
    while (rindicate==true){ // ind.expl : rejection sampling with exponential proposal. Use if lower >> mean
      a = (lower + std::sqrt(std::pow(lower,2.0)+4))/2.0;
      z = R::rexp(a) + lower;
      u = R::runif(0,1);
      if ((u <= std::exp(-std::pow(z-a,2.0)/2.0))&&(z <= upper)){
        rindicate = false;
        output = z;
      }
    }
  } else if ((upper<=0)&&(-lower>-upper+2.0*std::sqrt(std::exp(1.0))/(-upper + std::sqrt(std::pow(upper,2.0)+4))*std::exp((upper*2.0+upper*std::sqrt(std::pow(upper,2.0)+4))/4))){
    while (rindicate==true){ // ind.expu : rejection sampling with exponential proposal. Use if upper << mean.
      a = (-upper + std::sqrt(std::pow(upper,2.0)+4))/2.0;
      z = R::rexp(a)-upper;
      u = R::runif(0,1);
      if ((u<=std::exp(-std::pow(z-a,2.0)/2.0))&&(z<=-lower)){
        rindicate = false;
        output = -z;
      }
    }
  } else {
    while (rindicate==true){ // ind.u    : rejection sampling with uniform proposal. Use if bounds are narrow and central.
      z = R::runif(lower,upper);
      if (lower > 0){
        rho = std::exp((std::pow(lower, 2.0) - std::pow(z, 2.0))/2.0);
      } else {
        if (upper < 0){
          rho = std::exp((std::pow(upper,2.0)-std::pow(z,2.0))/2.0);
        } else {
          rho = std::exp(-(z*z)/2.0);
        }
      }
      u = R::runif(0,1);
      if (u<=rho){
        rindicate = false;
        output = z;
      }
    }
  }
  return(output*sd+mean);
}

// 2. dtruncnorm
// [[Rcpp::export]]
double dtruncnorm(double x, double mean, double sd, double lower, double upper){
  double output = 0.0, denom = 0.0, xtmp = 0.0;
  if ((x<lower)||(x>upper)){
    output = -arma::datum::inf;
  } else {
    denom  = R::pnorm(upper, mean, sd, 1, 0) - R::pnorm(lower, mean, sd, 1, 0); // 
    xtmp   = R::dnorm(x, mean, sd, 0);
    output = xtmp/denom;
  }
  return(output);
}
