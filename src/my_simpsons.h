#ifndef _carfima_MY_SIMPSONS_H
#define _carfima_MY_SIMPSONS_H
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include <cmath>
#include <complex>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;


double g1_reval(double x, arma::vec params);
double g1_ieval(double x, arma::vec params);
double g2_reval(double x, arma::vec params);
double g2_ieval(double x, arma::vec params);
  
double adaptiveSimpsonsAux(double (*f)(double, arma::vec), double a, double b, double epsilon,                 
                             double S, double fa, double fb, double fc, int bottom, arma::vec params);
double adaptiveSimpsonsAux(double (*f)(double, arma::vec), double a, double b, double epsilon,                 
                           double S, double fa, double fb, double fc, int bottom, arma::vec params);


double alt_g1real(double lower, double upper, double h_val, std::complex<double> d, double H);
double alt_g1imag(double lower, double upper, double h_val, std::complex<double> d, double H);
double alt_g2real(double lower, double upper, double h_val, std::complex<double> d, double H);
double alt_g2imag(double lower, double upper, double h_val, std::complex<double> d, double H);

#endif
