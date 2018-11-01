#ifndef _carfima_AUXILIARY_H
#define _carfima_AUXILIARY_H
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include <cmath>
#include <complex>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;
  
arma::cx_mat g1(arma::vec h, arma::cx_vec d, double H);
arma::cx_mat g2(arma::vec h, arma::cx_vec d, double H);

arma::mat B_mat(int p, arma::vec alpha);
arma::mat V_mat_sigma(int p, arma::vec alpha, double sigma);
arma::mat V_mat(int p, arma::vec alpha);

arma::cx_mat gamma_y_multiplication(double x, arma::cx_vec eigval, arma::cx_mat M2);
std::complex<double> cx_quadratic(arma::cx_vec xvec, arma::cx_mat xmat);
arma::cx_vec Gamma_Y(arma::vec time_lag_cov, int p, arma::mat A, double H, arma::cx_vec beta, arma::vec delta_p);
arma::cx_vec Gamma_Y_sigma(arma::vec time_lag_cov, int p, arma::mat A, double H, arma::cx_vec beta, arma::vec delta_p, double sigma);


arma::vec infvec(int n);
arma::mat infmat(int m, int n);
#endif
