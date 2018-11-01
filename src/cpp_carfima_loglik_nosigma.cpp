/*
 * MAIN FUNCTION FOR CARFIMA_LOGLIK_NOSIGMA
 *    - automatic segregation of parameters are already covered. 
 *    - I removed the 'fitted' boolean option; just return everything
 */

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "auxiliary.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List cpp_carfima_loglik_nosigma(arma::vec Yold, arma::vec time, int p, int q, arma::vec parameter){
  arma::vec alpha;
  arma::cx_vec beta;
  int ntime = time.n_elem;
  double H;
  if (q!=0L){ 
    alpha = parameter.subvec(0,(p-1));
    beta  = zeros<arma::cx_vec>(1+q+p-q-1);
    beta(0) = 1;
    for (int i=0;i<q;i++){
      beta(i+1) = parameter(p+i);
    }
    H = parameter(p+q);
  } else { 
    alpha = parameter.subvec(0,(p-1));
    beta  = zeros<arma::cx_vec>(p);
    beta(0) = 1;
    H = parameter(p);
  }
  arma::vec delta_p(p,fill::zeros);
  delta_p(p-1) = 1; 
  
  double mY = mean(Yold);
  arma::vec Y(Yold.n_elem, fill::zeros);
  for (int i=0;i<Yold.n_elem;i++){
    Y(i) = Yold(i)-mY;
  }
  
  arma::mat time_lag(time.n_elem, time.n_elem, fill::zeros);
  for (int i=0;i<time.n_elem;i++){
    for (int j=0;j<time.n_elem;j++){
      time_lag(i,j) = arma::datum::inf;
    }
  }
  for (int j=(time.n_elem-1);j>0;--j){
    for (int i=j;i>0;--i){
      time_lag(j,i-1) = time(j)-time(i-1);
    }
  } 
  arma::vec time_lag_cov = arma::sort(arma::vectorise(time_lag.elem(find_finite(time_lag)))); // ascending of non NAN elements in the matrix
  arma::mat A(p,p,fill::zeros);
  if (p!=1){
    for (int i=0;i<(p-1);i++){
      A(i,i+1) = 1;
    }
  }
  A.row(p-1) = alpha;
  
  arma::cx_vec Gamma_Y_cx = Gamma_Y(time_lag_cov, p, A, H, beta, delta_p);
  arma::vec Gamma_Y = arma::real(Gamma_Y_cx);

  return(Rcpp::List::create(Rcpp::Named("GammaY") = Gamma_Y,
                            Rcpp::Named("Y") = Y,
                            Rcpp::Named("time_lag_cov") = time_lag_cov));
}
