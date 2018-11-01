#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "loadcubature.h"

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _carfima_B_mat(SEXP, SEXP);
extern SEXP _carfima_cpp_carfima_bayes(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _carfima_cpp_carfima_loglik(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _carfima_cpp_carfima_loglik_nosigma(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _carfima_dtruncnorm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _carfima_g1(SEXP, SEXP, SEXP);
extern SEXP _carfima_g1_integration_imag(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _carfima_g1_integration_real(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _carfima_g2(SEXP, SEXP, SEXP);
extern SEXP _carfima_g2_integration_imag(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _carfima_g2_integration_real(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _carfima_Gamma_Y(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _carfima_Gamma_Y_sigma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _carfima_rtruncnorm(SEXP, SEXP, SEXP, SEXP);
extern SEXP _carfima_V_mat(SEXP, SEXP);
extern SEXP _carfima_V_mat_sigma(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_carfima_B_mat",                      (DL_FUNC) &_carfima_B_mat,                      2},
  {"_carfima_cpp_carfima_bayes",          (DL_FUNC) &_carfima_cpp_carfima_bayes,          9},
  {"_carfima_cpp_carfima_loglik",         (DL_FUNC) &_carfima_cpp_carfima_loglik,         5},
  {"_carfima_cpp_carfima_loglik_nosigma", (DL_FUNC) &_carfima_cpp_carfima_loglik_nosigma, 5},
  {"_carfima_dtruncnorm",                 (DL_FUNC) &_carfima_dtruncnorm,                 5},
  {"_carfima_g1",                         (DL_FUNC) &_carfima_g1,                         3},
  {"_carfima_g1_integration_imag",        (DL_FUNC) &_carfima_g1_integration_imag,        5},
  {"_carfima_g1_integration_real",        (DL_FUNC) &_carfima_g1_integration_real,        5},
  {"_carfima_g2",                         (DL_FUNC) &_carfima_g2,                         3},
  {"_carfima_g2_integration_imag",        (DL_FUNC) &_carfima_g2_integration_imag,        5},
  {"_carfima_g2_integration_real",        (DL_FUNC) &_carfima_g2_integration_real,        5},
  {"_carfima_Gamma_Y",                    (DL_FUNC) &_carfima_Gamma_Y,                    6},
  {"_carfima_Gamma_Y_sigma",              (DL_FUNC) &_carfima_Gamma_Y_sigma,              7},
  {"_carfima_rtruncnorm",                 (DL_FUNC) &_carfima_rtruncnorm,                 4},
  {"_carfima_V_mat",                      (DL_FUNC) &_carfima_V_mat,                      2},
  {"_carfima_V_mat_sigma",                (DL_FUNC) &_carfima_V_mat_sigma,                3},
  {NULL, NULL, 0}
};

void R_init_carfima(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  
  my_pcubature = (int (*) (unsigned, integrand, void *, unsigned, const double *, const double *,
                       size_t, double, double, error_norm, double *, double *)) R_GetCCallable("cubature","pcubature");
  my_hcubature = (int (*) (unsigned, integrand, void *, unsigned, const double *, const double *,
                       size_t, double, double, error_norm, double *, double *)) R_GetCCallable("cubature","hcubature");
}
