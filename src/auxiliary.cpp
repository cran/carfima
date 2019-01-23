// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <cmath>
#include <complex>
#include "loadcubature.h"
#include "cubature.h"
#include "auxiliary.h"
#include "my_simpsons.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

/*  
 *  (0) Intergration-Related Files
 *  (1) g1, g2
 *  (2) B.mat, V.mat.sigma, V.mat
 *  (3) Gamma.Y, Gamma.Y.sigma
 *  (4) infvec, infmat
 * 
 */
////////////////////////////////////////////////////////////////////////////////
// (0) integration-related files

int g1_integration_real_eval(unsigned dim, const double *u, void *parameters, unsigned fdim, double *fval){
  double *params = (double *) parameters;
  double h_val = params[0];
  double d_real = params[1];
  double d_imag = params[2];
  double H = params[3];
  
  fval[0] = (std::exp(d_real*(h_val- u[0])))*(std::cos(d_imag*(h_val-u[0])))*(std::pow(u[0],(2*H-2)));
  return(0);
}
int g1_integration_imag_eval(unsigned dim, const double *u, void *parameters, unsigned fdim, double *fval){
  double *params = (double *) parameters;
  double h_val = params[0];
  double d_real = params[1];
  double d_imag = params[2];
  double H = params[3];
  
  fval[0] = (std::exp(d_real*(h_val- u[0])))*(std::sin(d_imag*(h_val-u[0])))*(std::pow(u[0],(2*H-2)));
  return(0);
}

int g2_integration_real_eval(unsigned dim, const double *u, void *parameters, unsigned fdim, double *fval){
  double *params = (double *) parameters;
  double h_val = params[0];
  double d_real = params[1];
  double d_imag = params[2];
  double H = params[3];
  
  fval[0] = (std::exp(d_real*(u[0]-h_val)))*(std::cos(d_imag*(u[0]-h_val)))*(std::pow(u[0],(2*H-2)));
  return(0);
}
int g2_integration_imag_eval(unsigned dim, const double *u, void *parameters, unsigned fdim, double *fval){
  double *params = (double *) parameters;
  double h_val = params[0];
  double d_real = params[1];
  double d_imag = params[2];
  double H = params[3];
  
  fval[0] = (std::exp(d_real*(u[0]-h_val)))*(std::sin(d_imag*(u[0]-h_val)))*(std::pow(u[0],(2*H-2)));
  return(0);
}


// [[Rcpp::export]]
double g1_integration_real(double lower,double upper,double h_val,std::complex<double> d,double H){
  double params[4];
  params[0] = h_val;
  params[1] = d.real();
  params[2] = d.imag();
  params[3] = H;

  double val, err; // to be used from cubature call
  my_pcubature(1,g1_integration_real_eval,params,1,&lower,&upper,0,0.0001220703,0.0001220703,ERROR_INDIVIDUAL,&val,&err);
  return(val);
}
// [[Rcpp::export]]
double g1_integration_imag(double lower,double upper,double h_val,std::complex<double> d,double H){
  double params[4];
  params[0] = h_val;
  params[1] = d.real();
  params[2] = d.imag();
  params[3] = H;
  
  double val, err; // to be used from cubature call
  my_pcubature(1,g1_integration_imag_eval,params,1,&lower,&upper,0,0.00012207031,0.00012207030,ERROR_INDIVIDUAL,&val,&err);
  return(val);
}
// [[Rcpp::export]]
double g2_integration_real(double lower,double upper,double h_val,std::complex<double> d,double H){
  double params[4];
  params[0] = h_val;
  params[1] = d.real();
  params[2] = d.imag();
  params[3] = H;
  
  double val, err; // to be used from cubature call
  my_pcubature(1,g2_integration_real_eval,params,1,&lower,&upper,0,0.0001220703,0.0001220703,ERROR_INDIVIDUAL,&val,&err);
  return(val);
}
// [[Rcpp::export]]
double g2_integration_imag(double lower,double upper,double h_val,std::complex<double> d,double H){
  double params[4];
  params[0] = h_val;
  params[1] = d.real();
  params[2] = d.imag();
  params[3] = H;
  
  double val, err; // to be used from cubature call
  my_pcubature(1,g2_integration_imag_eval,params,1,&lower,&upper,0,0.0001220703,0.0001220703,ERROR_INDIVIDUAL,&val,&err);
  return(val);
}


double g1_integration_realh(double lower,double upper,double h_val,std::complex<double> d,double H){
  double params[4];
  params[0] = h_val;
  params[1] = d.real();
  params[2] = d.imag();
  params[3] = H;
  
  double val, err; // to be used from cubature call
  my_hcubature(1,g1_integration_real_eval,params,1,&lower,&upper,0,0.0001220703,0.0001220703,ERROR_INDIVIDUAL,&val,&err);
  return(val);
}
double g1_integration_imagh(double lower,double upper,double h_val,std::complex<double> d,double H){
  double params[4];
  params[0] = h_val;
  params[1] = d.real();
  params[2] = d.imag();
  params[3] = H;
  
  double val, err; // to be used from cubature call
  my_hcubature(1,g1_integration_imag_eval,params,1,&lower,&upper,0,0.00012207031,0.00012207030,ERROR_INDIVIDUAL,&val,&err);
  return(val);
}
double g2_integration_realh(double lower,double upper,double h_val,std::complex<double> d,double H){
  double params[4];
  params[0] = h_val;
  params[1] = d.real();
  params[2] = d.imag();
  params[3] = H;
  
  double val, err; // to be used from cubature call
  my_hcubature(1,g2_integration_real_eval,params,1,&lower,&upper,0,0.0001220703,0.0001220703,ERROR_INDIVIDUAL,&val,&err);
  return(val);
}
double g2_integration_imagh(double lower,double upper,double h_val,std::complex<double> d,double H){
  double params[4];
  params[0] = h_val;
  params[1] = d.real();
  params[2] = d.imag();
  params[3] = H;
  
  double val, err; // to be used from cubature call
  my_hcubature(1,g2_integration_imag_eval,params,1,&lower,&upper,0,0.0001220703,0.0001220703,ERROR_INDIVIDUAL,&val,&err);
  return(val);
}


////////////////////////////////////////////////////////////////////////////////
// (1) g1, g2
// [[Rcpp::export]]
arma::cx_mat g1(arma::vec h, arma::cx_vec d, double H){       // g1 function
  int nh = h.n_elem;
  int leng_h = h.n_elem;
  int leng_d = d.n_elem;
  arma::vec h0((nh+1), fill::zeros); // c(0, h) is replicated
  for (int i=1;i<(nh+1);i++){
    h0(i) = h(i-1);
  }

  // some variables to be declared
  double real_part_integration = 0;
  double imag_part_integration = 0;
  
  arma::cx_mat g1_res((leng_h+1),(leng_d),fill::zeros);
  std::complex<double> valuetoshow;
  for (int j=0; j<leng_d; j++){
    for (int i=0; i<leng_h; i++){
      if (i<=2){
        real_part_integration = g1_integration_realh(h0(i),h0(i+1),h0(i+1),d(j),H); 
        imag_part_integration = g1_integration_imagh(h0(i),h0(i+1),h0(i+1),d(j),H);  
      } else {
        real_part_integration = g1_integration_real(h0(i),h0(i+1),h0(i+1),d(j),H); 
        imag_part_integration = g1_integration_imag(h0(i),h0(i+1),h0(i+1),d(j),H);  
      }

      std::complex<double> tmp2add(real_part_integration, imag_part_integration);
      valuetoshow = std::exp( d(j)*(h0(i+1)-h0(i)))*g1_res(i,j);
      g1_res((i+1),j) = valuetoshow + tmp2add; 
    }
  }
  return(g1_res);
}
// [[Rcpp::export]]
arma::cx_mat g2(arma::vec h, arma::cx_vec d, double H){
  int leng_h = h.n_elem;
  arma::vec h0((leng_h + 1), fill::zeros); // c(0, h) is replicated
  for (int i=1;i<(leng_h+1);i++){
    h0(i) = h(i-1);
  }
  int leng_h0 = h0.n_elem;
  arma::vec h0_pseudo_max((leng_h0)+1);
  for (int i=0;i<leng_h0;i++){
    h0_pseudo_max(i) = h0(i);
  }
  h0_pseudo_max(leng_h0) = h0(leng_h0-1)*100.0;
  int leng_d = d.n_elem;
  
  
  // some variables to be declared
  double real_part_integration = 0;
  double imag_part_integration = 0;
  
  arma::cx_mat g2_res((leng_h0+1),(leng_d),fill::zeros);
  double glower, gupper, gh_value;
  std::complex<double> gd;
  int i;
  for (int j=0; j<leng_d; j++){
    for (int k=0;k<leng_h0;k++){
      i = leng_h0-k-1;
      glower   = h0_pseudo_max(i);
      gupper   = h0_pseudo_max(i+1);
      gh_value = h0_pseudo_max(i);
      gd       = d(j);
      
      if (i<=2){
        real_part_integration = g2_integration_realh(glower,gupper,gh_value, gd, H);
        imag_part_integration = g2_integration_imagh(glower,gupper,gh_value, gd, H);  
      } else {
        real_part_integration = g2_integration_real(glower,gupper,gh_value, gd, H);
        imag_part_integration = g2_integration_imag(glower,gupper,gh_value, gd, H);  
      }
      
      
      
      std::complex<double> tmp2add(real_part_integration, imag_part_integration);
      g2_res(i,j) = std::exp( d(j)*( h0_pseudo_max(i+1)-h0_pseudo_max(i)))*g2_res((i+1),j) + tmp2add;
    }
  }
  g2_res.shed_row(leng_h0);
  return(g2_res);
}

////////////////////////////////////////////////////////////////////////////////
// (2) B.mat, V.mat.sigma, V.mat
// [[Rcpp::export]]
arma::mat B_mat(int p, arma::vec alpha){
  arma::mat B(p,p,fill::zeros);
  
  for (int i=0;i<p;i++){
    for (int j=0;j<p;j++){
      if (((2*j-i+1)<=p)&&((2*j-i+1)>=1)){
        B(i,j) = std::pow(-1.0,j-i)*alpha((2*j)-i);
      } else if ((2*j-i+1)==(p+1)){
        B(i,j) = std::pow(-1.0,j-i-1.0);
      }
    }
  }
  return(B);
}
// [[Rcpp::export]]
arma::mat V_mat_sigma(int p, arma::vec alpha, double sigma){
  arma::mat B = B_mat(p,alpha);
  arma::colvec delta_p(p,fill::zeros);
  delta_p(p-1) = 1.0;
  arma::colvec V_diag = (-(sigma*sigma)/2)*arma::inv(B)*delta_p;
  
  arma::mat V(p,p,fill::zeros);
  for (int i=0;i<p;i++){
    for (int j=0;j<p;j++){
      if ((i+j)%2L == 0){
        V(i,j) = std::pow(-1.0, (static_cast<double>(i-j))/2.0)*V_diag(((i+j)/2));
      } else {
        V(i,j) = 0.0;
      }
    }
  }
  for (int u=0;u<p;u++){
    V(u,u) = V_diag(u);
  }
  return(V);
}
// [[Rcpp::export]]
arma::mat V_mat(int p, arma::vec alpha){
  arma::mat B = B_mat(p, alpha); 
  arma::colvec delta_p(p,fill::zeros);
  delta_p(p-1) = 1;
  arma::colvec V_diag = (-0.5)*arma::inv(B)*delta_p;
  
  arma::mat V(p,p,fill::zeros);
  for (int i=0;i<p;i++){
    for (int j=0;j<p;j++){
      if ((i+j)%2L == 0){
        V(i,j) = std::pow(-1.0, static_cast<double>(i-j)/2.0)*V_diag(((i+j)/2));
      }
    }
  }
  
  for (int u=0;u<p;u++){
    V(u,u) = V_diag(u);
  }
  return(V);
}


////////////////////////////////////////////////////////////////////////////////
// (3) Gamma.Y, Gamma.Y.sigma
arma::cx_mat gamma_y_multiplication(double x, arma::cx_vec eigval, arma::cx_mat M2){
  const int n = eigval.n_elem;
  arma::cx_vec multiplier(n,fill::zeros);
  for (int i=0;i<n;i++){
    multiplier(i) = std::exp(x*eigval(i));
  }
  arma::cx_mat output(n,n,fill::zeros);
  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++){
      output(i,j) = multiplier(i)*M2(i,j);
    }
  }
  return(output);
}
std::complex<double> cx_quadratic(arma::cx_vec xvec, arma::cx_mat xmat){
  const int n = xvec.n_elem;
  // arma::cx_colvec xcolvec = xvec;
  // arma::cx_rowvec xrowvec = arma::conv_to<cx_rowvec>::from(xcolvec); Rcpp::Rcout << "line cx_quadratic:331 complete" << std::endl;
  //   
  // std::complex<double> output = sum(xrowvec*xmat*xcolvec); Rcpp::Rcout << "line cx_quadratic:333 complete" << std::endl;
  // std::complex<double> output = sum(xvec.t()*xmat*xvec);
  std::complex<double> output = arma::as_scalar(xvec.t()*xmat*xvec); 
  return(output);
}
// [[Rcpp::export]]
arma::cx_vec Gamma_Y(arma::vec time_lag_cov, int p, arma::mat A, double H, arma::cx_vec beta, arma::vec delta_p){
  int tlc = time_lag_cov.n_elem; // length of time_lag_cov
  double num_s = NA_REAL;
  arma::cx_vec Gamma_Y_out((time_lag_cov.n_elem + 1), fill::zeros);
  for (int i=0;i<Gamma_Y_out.n_elem;i++){
    Gamma_Y_out(i) = num_s;
  }
  
  double C_H = H*(2*H - 1.0);          
  arma::rowvec alpha  = A.row(p-1);     
  arma::vec    alpha2 = arma::conv_to<vec>::from(alpha); 
  arma::cx_vec eig_val;
  arma::cx_mat eig_vec;
  eig_gen(eig_val, eig_vec, A);        
  arma::cx_mat eig_vec_inv = arma::inv((eig_vec)); 
  arma::vec time_lag_cov_0((tlc+1), fill::zeros);
  for (int i=1;i<(tlc+1);i++){
    time_lag_cov_0(i) = time_lag_cov(i-1);
  } 
  
  arma::mat V_ast = V_mat(p, alpha2); 
  
  arma::cx_mat g1_save = g1(time_lag_cov, eig_val, H); 
  arma::cx_mat g2_save = g2(time_lag_cov, eig_val, H); 
  
  arma::cx_cube M1(p,p,(tlc + 1), fill::zeros);
  arma::cx_cube M2(p,p,(tlc + 1), fill::zeros);
  for (int i=0;i<p;i++){
    for (int j=0;j<p;j++){
      for (int k=0;k<M1.n_slices;k++){
        M1(i,j,k) = num_s;
        M2(i,j,k) = num_s;
      }
    }
  }
  
  for (int i=0;i<(time_lag_cov.n_elem+1);i++){
    if (p < 2){
      M1.slice(i) = g1_save.row(i);
    } else {
      M1.slice(i) = arma::diagmat(g1_save.row(i));
    }
  }
  
  for (int i=0;i<(time_lag_cov.n_elem+1);i++){
    if (p < 2){
      M2.slice(i) = g2_save.row(i);
    } else {
      M2.slice(i) = arma::diagmat(g2_save.row(i));
    }
  }
  
  std::complex<double> gammaval(0,0);
  arma::cx_mat quadmat;
  for (int i=0;i<(time_lag_cov.n_elem+1);i++){
    quadmat = eig_vec*(M1.slice(i) + M2.slice(i) + gamma_y_multiplication(time_lag_cov_0(i), eig_val, M2.slice(0)))*eig_vec_inv*V_ast;
    // gammaval = cx_quadratic(beta, quadmat); 
    gammaval = arma::as_scalar(beta.t()*quadmat*beta);
    Gamma_Y_out(i) = (C_H * gammaval); 
  }
  return(Gamma_Y_out);
}
// [[Rcpp::export]]
arma::cx_vec Gamma_Y_sigma(arma::vec time_lag_cov, int p, arma::mat A, double H, arma::cx_vec beta, arma::vec delta_p, double sigma){
  int tlc = time_lag_cov.n_elem; // length of time_lag_cov
  double num_s = NA_REAL;
  arma::cx_vec Gamma_Y_out((time_lag_cov.n_elem + 1), fill::zeros);
  for (int i=0;i<Gamma_Y_out.n_elem;i++){
    Gamma_Y_out(i) = num_s;
  }
  double C_H = H*(2*H - 1.0);          
  arma::rowvec alpha  = A.row(p-1);    
  arma::vec    alpha2 = arma::conv_to<vec>::from(alpha);
  arma::vec eig_val;
  arma::mat eig_vec;
  eig_sym(eig_val, eig_vec, A);        
  arma::cx_mat eig_vec_inv = arma::conv_to<cx_mat>::from(arma::inv((eig_vec))); 
  arma::vec time_lag_cov_0((tlc+1), fill::zeros);
  for (int i=1;i<(tlc+1);i++){
    time_lag_cov_0(i) = time_lag_cov(i-1);
  } 
  
  arma::mat V_ast = V_mat_sigma(p, alpha2, sigma);
  
  arma::cx_vec cxeig_val = arma::conv_to<cx_vec>::from(eig_val);
  
  arma::cx_mat g1_save = g1(time_lag_cov, cxeig_val, H);
  arma::cx_mat g2_save = g2(time_lag_cov, cxeig_val, H);
  
  arma::cx_cube M1(p,p,(tlc + 1), fill::zeros);
  arma::cx_cube M2(p,p,(tlc + 1), fill::zeros);
  for (int i=0;i<p;i++){
    for (int j=0;j<p;j++){
      for (int k=0;k<M1.n_slices;k++){
        M1(i,j,k) = num_s;
        M2(i,j,k) = num_s;
      }
    }
  }
  for (int i=0;i<(time_lag_cov.n_elem+1);i++){
    if (p < 2){
      M1.slice(i) = sum(g1_save(i));
    } else {
      M1.slice(i) = arma::diagmat(g1_save.row(i));
    }
  }
  for (int i=0;i<(time_lag_cov.n_elem+1);i++){
    if (p < 2){
      M2.slice(i) = sum(g2_save(i));
    } else {
      M2.slice(i) = arma::diagmat(g2_save.row(i));
    }
  }
  

  
  std::complex<double> gammaval(0,0);
  arma::cx_mat quadmat;
  for (int i=0;i<(time_lag_cov.n_elem+1);i++){
    quadmat = eig_vec*(M1.slice(i) + M2.slice(i) + gamma_y_multiplication(time_lag_cov_0(i), cxeig_val, M2.slice(0)))*eig_vec_inv*V_ast;
    gammaval = arma::as_scalar(beta.t()*quadmat*beta);
    // cx_quadratic(beta, quadmat);
    Gamma_Y_out(i) = (C_H * gammaval);
  }
  return(Gamma_Y_out);
}
////////////////////////////////////////////////////////////////////////////////
// (4) infvec, infmat
arma::vec infvec(int n){
  arma::vec output(n,fill::zeros);
  for (int i=0;i<n;i++){
    output(i) = arma::datum::inf;
  }
  return(output);
}
arma::mat infmat(int m, int n){
  arma::mat output(m,n,fill::zeros);
  for (int i=0;i<m;i++){
    for (int j=0;j<m;j++){
      output(i,j) = arma::datum::inf;
    }
  }
  return(output);
}
