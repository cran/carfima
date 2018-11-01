#include "cubature.h"

int (*my_pcubature)(unsigned, integrand, void *, unsigned, const double *, const double *,
     size_t, double, double, error_norm, double *, double *);
int (*my_hcubature)(unsigned, integrand, void *, unsigned, const double *, const double *,
     size_t, double, double, error_norm, double *, double *);
