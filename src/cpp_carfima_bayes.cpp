#include <RcppArmadillo.h>
#include "auxiliary.h"
#include "my_truncnorm.h"

// [[Rcpp::depends(RcppArmadillo)]]

// use of truncnorm
// double rtruncnorm(double mean, double sd, double lbd, double rbd)
// double dtruncnorm(double x, double mean, double sd, double lower, double upper)

// extra function : dinvgamma
// carfima.loglik.internal <- function(Y, time, ar.p, ma.q, parameter) will be used

double my_dinvgamma(double x, double shape, double scale, bool log){
  double rate = scale;
  double log_f = R::dgamma(1.0/x, shape, rate, 1) - 2.0*static_cast<double>(std::log(static_cast<float>(x)));
  if (log){
    return(log_f);
  } else {
    return(static_cast<double>(std::exp(static_cast<float>(log_f))));
  }
}
double my_log(double x){
  return(static_cast<double>(std::log(static_cast<float>(x))));
}
double my_exp(double x){
  return(static_cast<double>(std::exp(static_cast<float>(x))));
}
double my_sqrt(double x){
  return(static_cast<double>(std::sqrt(static_cast<float>(x))));
}



// [[Rcpp::export]]
Rcpp::List cpp_carfima_bayes(Function loglik, arma::vec Y, arma::vec time, 
                             arma::vec param_ini, arma::vec param_scale,  int ar_p, int ma_q,
                             int n_warm, int n_sample){
  
  /*
   * PREPROCESSING : reproduce initial parameters 
   */
  int p = ar_p;
  int q = ma_q;
  int n_total = n_warm + n_sample;
  int dimen   = param_ini.n_elem;
  
  arma::mat accept_out(n_total, dimen, fill::zeros);
  arma::mat param_out(n_total, dimen, fill::zeros);
  arma::vec param_t = param_ini;
  
  double prev_log_den = as<double>(loglik(Y, time, p, q, param_t));
  double meanapt = 0.0;
  double sqrti00 = 0.0; 
  
  /*
   * MAIN ITERATION
   */
  // Part 0. Initialization
  double scale_adj;
  double param_p, prop_log_den;
  double l_metro, l_hastings, l_threshold;
  arma::vec param_temp, param_p_vec;
  
  for (int i=0;i<n_total;i++){
    // Part 1. Iterate over 'p'
    for (int j=0;j<p;j++){
      param_p = rtruncnorm(param_t(j), param_scale(j), -0.99, -0.01);
      param_temp = param_t;
      param_temp(j) = param_p;
      param_p_vec = param_temp;
      
      prop_log_den = as<double>(loglik(Y,time,p,q,param_p_vec));
      
      l_metro = prop_log_den - prev_log_den;
      l_hastings = my_log(dtruncnorm(param_t(j),param_p,param_scale(j),-0.99,-0.01))-my_log(dtruncnorm(param_p,param_t(j),param_scale(j),-0.99,-0.01));
      
      l_threshold = -R::rexp(1.0);
      if ((l_metro + l_hastings) > l_threshold){
        param_t = param_p_vec;
        prev_log_den = prop_log_den;
        accept_out(i,j) = 1.0;
      }
      
      meanapt = 0.0;
      sqrti00 = (1.0/(my_sqrt(static_cast<double>(i+1)/100.0)));
      arma::vec vecj;
      if (i>0){ // line I added 
        if ((i+1)%100==0L){
          vecj    = accept_out.col(j);
          meanapt = arma::mean(vecj.subvec(i-99,i));
          
          if (meanapt > 0.3){
            if (0.1 < sqrti00){
              scale_adj = my_exp(0.1);
            } else {
              scale_adj = my_exp(sqrti00);
            }
          } else if (meanapt < 0.3){
            if (0.1 < sqrti00){
              scale_adj = my_exp(-0.1);
            } else {
              scale_adj = my_exp(-sqrti00);
            }
          } else {
            scale_adj = 1.0;
          }
          param_scale(j) = param_scale(j)*scale_adj;
        }
      }
    }
    
    // Part 2.
    if (ma_q!=0L){
      for (int j=0;j<q;j++){
        // ma parameter update
        param_p = R::rnorm(param_t(p+j), param_scale(p+j));
        param_temp = param_t;
        param_temp(p+j) = param_p;
        param_p_vec = param_temp;
        
        prop_log_den = as<double>(loglik(Y, time, p, q, param_p_vec));
        l_metro = prop_log_den - prev_log_den + (R::dnorm(param_p, 0, 10000, true) - R::dnorm(param_t(p+j), 0, 10000, true));
        
        if (l_metro > (-R::rexp(1.0))){
          param_t = param_p_vec;
          prev_log_den = prop_log_den;
          accept_out(i,p+j) = 1.0;
        }
        
        meanapt = 0.0;
        sqrti00 = (1.0/(my_sqrt(static_cast<double>(i+1)/100.0)));
        arma::vec vecj;
        if (i>0){ // line I added 
          if ((i+1)%100==0L){
            vecj    = accept_out.col(p+j);
            meanapt = arma::mean(vecj.subvec(i-99,i));
            
            if (meanapt > 0.3){
              if (0.1 < sqrti00){
                scale_adj = my_exp(0.1);
              } else {
                scale_adj = my_exp(sqrti00);
              }
            } else if (meanapt < 0.3){
              if (0.1 < sqrti00){
                scale_adj = my_exp(-0.1);
              } else {
                scale_adj = my_exp(-sqrti00);
              }
            } else {
              scale_adj = 1.0;
            }
            param_scale(p+j) = param_scale(p+j)*scale_adj;
          }
        }
      }
    }
    
    // Part 3. /////////////////////////// param.t(p+q+1).. hmm..
    param_p = rtruncnorm(param_t(p+q), param_scale(p+q), 0.51, 0.99);
    param_temp = param_t;
    param_temp(p+q) = param_p;
    param_p_vec = param_temp;
    
    prop_log_den = as<double>(loglik(Y, time, p, q, param_p_vec));
    l_metro = prop_log_den - prev_log_den;
    l_hastings = my_log(dtruncnorm(param_t(p+q),param_p,param_scale(p+q),0.51,0.99))-my_log(dtruncnorm(param_p,param_t(p+q),param_scale(p+q),0.51,0.99));
    
    if ((l_metro + l_hastings) > -(R::rexp(1.0))){
      param_t = param_p_vec;
      prev_log_den = prop_log_den;
      accept_out(i,p+q) = 1.0;
    }
    // Part 4.
    meanapt = 0.0;
    sqrti00 = (1.0/(my_sqrt(static_cast<double>(i+1)/100.0)));
    arma::vec vecj;
    if (i>0){ // line I added 
      if ((i+1)%100==0L){
        vecj    = accept_out.col(p+q);
        meanapt = arma::mean(vecj.subvec(i-99,i));
        
        if (meanapt > 0.3){
          if (0.1 < sqrti00){
            scale_adj = my_exp(0.1);
          } else {
            scale_adj = my_exp(sqrti00);
          }
        } else if (meanapt < 0.3){
          if (0.1 < sqrti00){
            scale_adj = my_exp(-0.1);
          } else {
            scale_adj = my_exp(-sqrti00);
          }
        } else {
          scale_adj = 1.0;
        }
        param_scale(p+q) = param_scale(p+q)*scale_adj;
      }
    }
    
    // Part 5. sigma update
    param_p = my_exp((R::rnorm((2.0*my_log(param_t(p+q+1))), param_scale(p+q+1))));
    param_temp = param_t;
    param_temp(p+q+1) = my_sqrt(param_p);
    param_p_vec = param_temp;
    
    prop_log_den = as<double>(loglik(Y, time, p, q, param_p_vec));
    l_metro    = prop_log_den + my_dinvgamma(param_p, 2.01, 1000, true) - prev_log_den - my_dinvgamma(std::pow(param_t(p+q+1), 2.0), 2.01, 1000, true);
    l_hastings = my_log(param_p) - 2.0*my_log(param_t(p+q+1));

    // Part 6. accept-reject
    if ((l_metro + l_hastings) > -R::rexp(1.0)){
      param_t = param_p_vec;
      prev_log_den = prop_log_den;
      accept_out(i,p+q+1) = 1.0;
    }
    
    meanapt = 0.0;
    sqrti00 = (1.0/(my_sqrt(static_cast<double>(i+1)/100.0)));
    vecj.reset();
    if (i>0){ // line I added 
      if ((i+1)%100==0L){
        vecj    = accept_out.col(p+q+1);
        meanapt = arma::mean(vecj.subvec(i-99,i));
        
        if (meanapt > 0.3){
          if (0.1 < sqrti00){
            scale_adj = my_exp(0.1);
          } else {
            scale_adj = my_exp(sqrti00);
          }
        } else if (meanapt < 0.3){
          if (0.1 < sqrti00){
            scale_adj = my_exp(-0.1);
          } else {
            scale_adj = my_exp(-sqrti00);
          }
        } else {
          scale_adj = 1.0;
        }
        param_scale(p+q+1) = param_scale(p+q+1)*scale_adj;
      }
    }
    param_out.row(i) = param_t.t();
  }
  // fake return
  Rcpp::List out;
  out["pre_param"] = param_out;
  out["pre_accept"] = accept_out;
  return(out);
}

