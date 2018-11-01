/*
 *  www2.math.umd.edu/~mariakc/teaching-2/adaptivesimpson.c
 */
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <cmath>
#include <complex>
#include "my_simpsons.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

double g1_reval(double x, arma::vec params){
  double h_val = params(0);
  double d_real = params(1);
  double d_imag = params(2);
  double H = params(3);
  
  return(std::exp(d_real*(h_val- x)))*(std::cos(d_imag*(h_val-x)))*(std::pow(x,(2*H-2)));
}
double g1_ieval(double x, arma::vec params){
  double h_val = params(0);
  double d_real = params(1);
  double d_imag = params(2);
  double H = params(3);
  
  return(std::exp(d_real*(h_val- x)))*(std::sin(d_imag*(h_val-x)))*(std::pow(x,(2*H-2)));
}
double g2_reval(double x, arma::vec params){
  double h_val = params(0);
  double d_real = params(1);
  double d_imag = params(2);
  double H = params(3);
  
  return(std::exp(d_real*(x-h_val)))*(std::cos(d_imag*(x-h_val)))*(std::pow(x,(2*H-2)));
}
double g2_ieval(double x, arma::vec params){
  double h_val = params(0);
  double d_real = params(1);
  double d_imag = params(2);
  double H = params(3);
  
  return(std::exp(d_real*(x-h_val)))*(std::sin(d_imag*(x-h_val)))*(std::pow(x,(2*H-2)));
}


double adaptiveSimpsonsAux(double (*f)(double, arma::vec), double a, double b, double epsilon,                 
                           double S, double fa, double fb, double fc, int bottom, arma::vec params) {                 
  double c = (a + b)/2, h = b - a;                                                                  
  double d = (a + c)/2, e = (c + b)/2;                                                              
  double fd = f(d, params), fe = f(e, params);                                                                      
  double Sleft = (h/12)*(fa + 4*fd + fc);                                                           
  double Sright = (h/12)*(fc + 4*fe + fb);                                                          
  double S2 = Sleft + Sright;                                                                       
  if (bottom <= 0 || std::abs(S2 - S) <= 15*epsilon) {
    if( bottom <=0 ) 
    return S2 + (S2 - S)/15;      
  }
  return adaptiveSimpsonsAux(f, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1, params) +                    
    adaptiveSimpsonsAux(f, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1, params);                     
}  

double adaptiveSimpsons(double (*f)(double, arma::vec),   // ptr to function
                        double a, double b,  // interval [a,b]
                        double epsilon,  // error tolerance
                        int maxRecursionDepth, arma::vec params){
  double c = (a + b)/2, h = b - a;                                                                  
  double fa = f(a, params), fb = f(b, params), fc = f(c, params);                                                           
  double S = (h/6)*(fa + 4*fc + fb); 
  return adaptiveSimpsonsAux(f, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth, params); 
}


double alt_g1real(double lower, double upper, double h_val, std::complex<double> d, double H){
  arma::vec params(4,fill::zeros);
  params(0) = h_val;
  params(1) = d.real();
  params(2) = d.imag();
  params(3) = H;
  
  return(adaptiveSimpsons(g1_reval, lower, upper, 1.0e-6, 40, params));
}
double alt_g1imag(double lower, double upper, double h_val, std::complex<double> d, double H){
  arma::vec params(4,fill::zeros);
  params(0) = h_val;
  params(1) = d.real();
  params(2) = d.imag();
  params(3) = H;
  
  return(adaptiveSimpsons(g1_ieval, lower, upper, 1.0e-6, 40, params));
}
double alt_g2real(double lower, double upper, double h_val, std::complex<double> d, double H){
  arma::vec params(4,fill::zeros);
  params(0) = h_val;
  params(1) = d.real();
  params(2) = d.imag();
  params(3) = H;
  
  return(adaptiveSimpsons(g2_reval, lower, upper, 1.0e-6, 40, params));
}
double alt_g2imag(double lower, double upper, double h_val, std::complex<double> d, double H){
  arma::vec params(4,fill::zeros);
  params(0) = h_val;
  params(1) = d.real();
  params(2) = d.imag();
  params(3) = H;
  
  return(adaptiveSimpsons(g2_ieval, lower, upper, 1.0e-6, 40, params));
}
